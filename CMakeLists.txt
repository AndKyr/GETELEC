# Set the minimum version of CMake required
cmake_minimum_required(VERSION 3.10)

# Define the project name and version
project(GETELEC VERSION 3.0 LANGUAGES CXX)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)


# Add an option to set the build type (default is Release if not specified)
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type (default is Release)" FORCE)
endif()


# Include directories
include_directories(${PROJECT_SOURCE_DIR}/include)

# Find GSL library
find_package(GSL REQUIRED)

# Find TBB library
find_package(TBB REQUIRED)


# Add the source files
file(GLOB_RECURSE SOURCES ${PROJECT_SOURCE_DIR}/src/*.cpp)

# Define a library target for your project
add_library(getelec_lib SHARED ${SOURCES})

# Link the GSL and math library
target_link_libraries(getelec_lib TBB::tbb GSL::gsl GSL::gslcblas m)

# Set the library's output name to "getelec"
set_target_properties(getelec_lib PROPERTIES OUTPUT_NAME getelec)

# Create the executable
add_executable(getelec ${PROJECT_SOURCE_DIR}/src/main.cpp)

# Link the GSL and math library
target_link_libraries(getelec getelec_lib)

# Add debug flags if building in Debug mode
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_definitions(getelec_lib PRIVATE DEBUG_MODE)
    target_compile_options(getelec_lib PRIVATE -g -O0 -DDEBUG)
    target_compile_definitions(getelec PRIVATE DEBUG_MODE)
    target_compile_options(getelec PRIVATE -g -O0 -DDEBUG)
endif()


# Add release-specific flags if building in Release mode
if (CMAKE_BUILD_TYPE STREQUAL "Release")
    target_compile_options(getelec_lib PRIVATE -O3 -DRELEASE)
    target_compile_options(getelec PRIVATE -O3 -DRELEASE)
endif()


#enable the tests
enable_testing()

# Include the tests directory
add_subdirectory(tests)

