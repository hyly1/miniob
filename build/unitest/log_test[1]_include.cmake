if(EXISTS "/Users/wangyunlai.wyl/Documents/project/miniob-ce/build/unitest/log_test[1]_tests.cmake")
  include("/Users/wangyunlai.wyl/Documents/project/miniob-ce/build/unitest/log_test[1]_tests.cmake")
else()
  add_test(log_test_NOT_BUILT log_test_NOT_BUILT)
endif()