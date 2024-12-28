# tests/CMakeLists.txt

# Add test executable
add_executable(tests test_cases.cpp)

# Link Google Test libraries
target_link_libraries(tests gtest gtest_main)

# Link your project libraries
target_include_directories(tests PRIVATE ${CMAKE_SOURCE_DIR}/include)

# Link Google Test libraries
target_link_libraries(tests gtest gtest_main getelec_lib GSL::gsl GSL::gslcblas m)

# target_link_libraries(tests PRIVATE getelec_lib)


# Add test command
add_test(NAME run_tests COMMAND tests)
