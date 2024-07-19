/* Copyright (c) 2021 OceanBase and/or its affiliates. All rights reserved.
miniob is licensed under Mulan PSL v2.
You can use this software according to the terms and conditions of the Mulan PSL v2.
You may obtain a copy of Mulan PSL v2 at:
         http://license.coscl.org.cn/MulanPSL2
THIS SOFTWARE IS PROVIDED ON AN "AS IS" BASIS, WITHOUT WARRANTIES OF ANY KIND,
EITHER EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO NON-INFRINGEMENT,
MERCHANTABILITY OR FIT FOR A PARTICULAR PURPOSE.
See the Mulan PSL v2 for more details. */

#include "sql/expr/aggregate_hash_table.h"
#include <iostream>

// ----------------------------------StandardAggregateHashTable------------------

RC StandardAggregateHashTable::add_chunk(Chunk &groups_chunk, Chunk &aggrs_chunk)
{
  if (groups_chunk.rows() != aggrs_chunk.rows()) {
    LOG_WARN("group_chunk and aggr_chunk rows must be equal.");
    return RC::INVALID_ARGUMENT;
  }
  for (int row = 0; row < groups_chunk.rows(); ++row) {
    std::vector<Value> group_values;
    for (int col = 0; col < groups_chunk.column_num(); ++col) {
      group_values.push_back(groups_chunk.column(col).get_value(row));
    }
    auto it = aggr_values_.find(group_values);
    if (it == aggr_values_.end()) {
      std::vector<Value> aggr_values;
      for (int col = 0; col < aggrs_chunk.column_num(); ++col) {
        aggr_values.push_back(aggrs_chunk.column(col).get_value(row));
      }
      aggr_values_[group_values] = aggr_values;
    } else {
      for (int col = 0; col < aggrs_chunk.column_num(); ++col) {
        auto &aggr_value = it->second[col];
        auto  cell_value = aggrs_chunk.column(col).get_value(row);
        switch (aggr_types_[col]) {
          case AggregateExpr::Type::SUM:
            if (aggr_value.attr_type() == AttrType::INTS && cell_value.attr_type() == AttrType::INTS) {
              aggr_value.set_int(aggr_value.get_int() + cell_value.get_int());
              LOG_WARN("value: %d", aggr_value.get_int());
            } else if (aggr_value.attr_type() == AttrType::FLOATS && cell_value.attr_type() == AttrType::FLOATS) {
              aggr_value.set_float(aggr_value.get_float() + cell_value.get_float());
            } else {
              LOG_WARN("unsupported type for SUM operation");
              return RC::INVALID_ARGUMENT;
            }
            break;
          case AggregateExpr::Type::MAX:
            if (aggr_value.compare(cell_value) < 0) {
              aggr_value.set_value(cell_value);
            }
            break;
          case AggregateExpr::Type::MIN:
            if (aggr_value.compare(cell_value) > 0) {
              aggr_value.set_value(cell_value);
            }
            break;
          case AggregateExpr::Type::COUNT: aggr_value.set_int(aggr_value.get_int() + 1); break;
          default: LOG_WARN("unsupported aggregate type"); return RC::INVALID_ARGUMENT;
        }
      }
    }
  }
  return RC::SUCCESS;
}

void StandardAggregateHashTable::Scanner::open_scan()
{
  it_  = static_cast<StandardAggregateHashTable *>(hash_table_)->begin();
  end_ = static_cast<StandardAggregateHashTable *>(hash_table_)->end();
}

RC StandardAggregateHashTable::Scanner::next(Chunk &output_chunk)
{
  if (it_ == end_) {
    return RC::RECORD_EOF;
  }
  while (it_ != end_ && output_chunk.rows() <= output_chunk.capacity()) {
    auto &group_by_values = it_->first;
    auto &aggrs           = it_->second;
    for (int i = 0; i < output_chunk.column_num(); i++) {
      auto col_idx = output_chunk.column_ids(i);
      if (col_idx >= static_cast<int>(group_by_values.size())) {
        output_chunk.column(i).append_one((char *)aggrs[col_idx - group_by_values.size()].data());
      } else {
        output_chunk.column(i).append_one((char *)group_by_values[col_idx].data());
      }
    }
    it_++;
  }
  if (it_ == end_) {
    return RC::SUCCESS;
  }

  return RC::SUCCESS;
}

size_t StandardAggregateHashTable::VectorHash::operator()(const vector<Value> &vec) const
{
  size_t hash = 0;
  for (const auto &elem : vec) {
    hash ^= std::hash<string>()(elem.to_string());
  }
  return hash;
}

bool StandardAggregateHashTable::VectorEqual::operator()(const vector<Value> &lhs, const vector<Value> &rhs) const
{
  if (lhs.size() != rhs.size()) {
    return false;
  }
  for (size_t i = 0; i < lhs.size(); ++i) {
    if (rhs[i].compare(lhs[i]) != 0) {
      return false;
    }
  }
  return true;
}

// ----------------------------------LinearProbingAggregateHashTable------------------
#ifdef USE_SIMD
template <typename V>
RC LinearProbingAggregateHashTable<V>::add_chunk(Chunk &group_chunk, Chunk &aggr_chunk)
{
  if (group_chunk.column_num() != 1 || aggr_chunk.column_num() != 1) {
    LOG_WARN("group_chunk and aggr_chunk size must be 1.");
    return RC::INVALID_ARGUMENT;
  }
  if (group_chunk.rows() != aggr_chunk.rows()) {
    LOG_WARN("group_chunk and aggr _chunk rows must be equal.");
    return RC::INVALID_ARGUMENT;
  }
  add_batch((int *)group_chunk.column(0).data(), (V *)aggr_chunk.column(0).data(), group_chunk.rows());
  return RC::SUCCESS;
}

template <typename V>
void LinearProbingAggregateHashTable<V>::Scanner::open_scan()
{
  capacity_   = static_cast<LinearProbingAggregateHashTable *>(hash_table_)->capacity();
  size_       = static_cast<LinearProbingAggregateHashTable *>(hash_table_)->size();
  scan_pos_   = 0;
  scan_count_ = 0;
}

template <typename V>
RC LinearProbingAggregateHashTable<V>::Scanner::next(Chunk &output_chunk)
{
  if (scan_pos_ >= capacity_ || scan_count_ >= size_) {
    return RC::RECORD_EOF;
  }
  auto linear_probing_hash_table = static_cast<LinearProbingAggregateHashTable *>(hash_table_);
  while (scan_pos_ < capacity_ && scan_count_ < size_ && output_chunk.rows() <= output_chunk.capacity()) {
    int key;
    V   value;
    RC  rc = linear_probing_hash_table->iter_get(scan_pos_, key, value);
    if (rc == RC::SUCCESS) {
      output_chunk.column(0).append_one((char *)&key);
      output_chunk.column(1).append_one((char *)&value);
      scan_count_++;
    }
    scan_pos_++;
  }
  return RC::SUCCESS;
}

template <typename V>
void LinearProbingAggregateHashTable<V>::Scanner::close_scan()
{
  capacity_   = -1;
  size_       = -1;
  scan_pos_   = -1;
  scan_count_ = 0;
}

template <typename V>
RC LinearProbingAggregateHashTable<V>::get(int key, V &value)
{
  RC  rc          = RC::SUCCESS;
  int index       = (key % capacity_ + capacity_) % capacity_;
  int iterate_cnt = 0;
  while (true) {
    if (keys_[index] == EMPTY_KEY) {
      rc = RC::NOT_EXIST;
      break;
    } else if (keys_[index] == key) {
      value = values_[index];
      break;
    } else {
      index += 1;
      index %= capacity_;
      iterate_cnt++;
      if (iterate_cnt > capacity_) {
        rc = RC::NOT_EXIST;
        break;
      }
    }
  }
  return rc;
}

template <typename V>
RC LinearProbingAggregateHashTable<V>::iter_get(int pos, int &key, V &value)
{
  RC rc = RC::SUCCESS;
  if (keys_[pos] == LinearProbingAggregateHashTable<V>::EMPTY_KEY) {
    rc = RC::NOT_EXIST;
  } else {
    key   = keys_[pos];
    value = values_[pos];
  }
  return rc;
}

template <typename V>
void LinearProbingAggregateHashTable<V>::aggregate(V *value, V value_to_aggregate)
{
  if (aggregate_type_ == AggregateExpr::Type::SUM) {
    *value += value_to_aggregate;
  } else {
    ASSERT(false, "unsupported aggregate type");
  }
}

template <typename V>
void LinearProbingAggregateHashTable<V>::resize()
{
  capacity_ *= 2;
  std::vector<int> new_keys(capacity_);
  std::vector<V>   new_values(capacity_);

  for (size_t i = 0; i < keys_.size(); i++) {
    auto &key   = keys_[i];
    auto &value = values_[i];
    if (key != EMPTY_KEY) {
      int index = (key % capacity_ + capacity_) % capacity_;
      while (new_keys[index] != EMPTY_KEY) {
        index = (index + 1) % capacity_;
      }
      new_keys[index]   = key;
      new_values[index] = value;
    }
  }

  keys_   = std::move(new_keys);
  values_ = std::move(new_values);
}

template <typename V>
void LinearProbingAggregateHashTable<V>::resize_if_need()
{
  if (size_ >= capacity_ / 2) {
    resize();
  }
}

// 一个函数用于将 __m256i 的内容输出到 std::ostream
void print_m256i(__m256i var)
{
  // 将 __m256i 转换为 8 个 32 位整数
  alignas(32) int32_t values[8];
  _mm256_store_si256(reinterpret_cast<__m256i *>(values), var);

  // 输出这些整数
  std::cout << "{ ";
  for (int i = 0; i < 8; ++i) {
    std::cout << values[i];
    if (i != 7) {
      std::cout << ", ";
    }
  }
  std::cout << " }" << std::endl;
}

template <typename V>
void LinearProbingAggregateHashTable<V>::add_batch(int *input_keys, V *input_values, int len)
{
  resize_if_need();
  __m256i inv = _mm256_set1_epi32(-1);
  __m256i off = _mm256_setzero_si256();
  __m256i hash;
  int     i          = 0;
  int     update_num = 0;
  int     keys[SIMD_WIDTH];
  V       values[SIMD_WIDTH];
  for (; i + SIMD_WIDTH <= len;) {
    // 1. Selective load

    selective_load(input_keys, i, keys, inv);
    selective_load(input_values, i, values, inv);

    // 2. Update i
    update_num = __builtin_popcountll(_mm256_movemask_ps(_mm256_castsi256_ps(inv)));
    i += update_num;

    // 3. Calculate hash
    hash =
        _mm256_and_si256(_mm256_set1_epi32(capacity_ - 1), _mm256_add_epi32(_mm256_loadu_si256((__m256i *)keys), off));

    // 4. Update aggregation results
    for (int j = 0; j < SIMD_WIDTH; j++) {
      int index = mm256_extract_epi32_var_indx(hash, j);
      if (keys_[index] == keys[j]) {
        aggregate(&values_[index], values[j]);
      } else if (keys_[index] == EMPTY_KEY) {
        keys_[index]   = keys[j];
        values_[index] = values[j];
        size_++;
      }
    }

    // 5. Gather keys from hash table
    __m256i table_key = _mm256_i32gather_epi32(keys_.data(), hash, 4);

    // 6. Update inv and off
    __m256i key_match = _mm256_cmpeq_epi32(table_key, _mm256_loadu_si256((__m256i *)keys));
    inv               = _mm256_blendv_epi8(_mm256_setzero_si256(), _mm256_set1_epi32(-1), key_match);
    // off               = _mm256_add_epi32(off, _mm256_andnot_si256(key_match, _mm256_set1_epi32(1)));
    __m256i add_one = _mm256_andnot_si256(key_match, _mm256_set1_epi32(1));
    __m256i new_off = _mm256_add_epi32(off, add_one);
    off             = _mm256_blendv_epi8(new_off, _mm256_setzero_si256(), key_match);

    // // Reset off for the next batch
    // if (_mm256_testz_si256(inv, inv)) {
    //   off = _mm256_setzero_si256();
    // }
  }
  // 处理剩余冲突
  alignas(32) int inv_array[SIMD_WIDTH];
  _mm256_store_si256((__m256i *)inv_array, inv);
  for (int j = 0; j < SIMD_WIDTH; j++) {
    if (inv_array[j] != -1) {
      int key   = keys[j];
      V   value = values[j];
      int index = key % (capacity_);
      while (true) {
        if (keys_[index] == EMPTY_KEY) {
          keys_[index]   = key;
          values_[index] = value;
          size_++;
          break;
        } else if (keys_[index] == key) {
          aggregate(&values_[index], value);
          break;
        }
        index = (index + 1) % capacity_;
      }
    }
  }
  // 7. Process remaining keys with scalar operations
  for (; i < len; i++) {
    int key   = input_keys[i];
    V   value = input_values[i];
    int index = key % (capacity_);
    while (true) {
      if (keys_[index] == EMPTY_KEY) {
        keys_[index]   = key;
        values_[index] = value;
        size_++;
        break;
      } else if (keys_[index] == key) {
        aggregate(&values_[index], value);
        break;
      }
      index = (index + 1) % capacity_;
    }
  }
}

template <typename V>
const int LinearProbingAggregateHashTable<V>::EMPTY_KEY = 0xffffffff;
template <typename V>
const int LinearProbingAggregateHashTable<V>::DEFAULT_CAPACITY = 16384;

template class LinearProbingAggregateHashTable<int>;
template class LinearProbingAggregateHashTable<float>;
#endif