# Set the minimum version of CMake required
cmake_minimum_required(VERSION 3.10)

# Define the project name and version
project(GETELEC VERSION 3.0 LANGUAGES CXX)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)


# Add an option to set the build type (default is Release if not specified)
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type (default is Release)" FORCE)
endif()


# Include directories
include_directories(${PROJECT_SOURCE_DIR}/include)

# Find GSL library
find_package(GSL REQUIRED)

# Add the source files
file(GLOB_RECURSE SOURCES ${PROJECT_SOURCE_DIR}/src/*.cpp)

# Create the executable
add_executable(getelec ${SOURCES})

# Link the GSL and math library
target_link_libraries(getelec GSL::gsl GSL::gslcblas m)

# Add debug flags if building in Debug mode
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_definitions(getelec PRIVATE DEBUG_MODE)
    target_compile_options(getelec PRIVATE -g -O0 -DDEBUG)
endif()


# Add release-specific flags if building in Release mode
if (CMAKE_BUILD_TYPE STREQUAL "Release")
    target_compile_options(getelec PRIVATE -O3 -DRELEASE)
endif()
