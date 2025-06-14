# Set the minimum version of CMake required
cmake_minimum_required(VERSION 3.10)

# Define the project name and version
project(GETELEC VERSION 3.0 LANGUAGES CXX C)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

option(ENABLE_JAVA_WRAPPER "Enable Java JNI Wrapper" OFF)


# Add an option to set the build type (default is Release if not specified)
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type (default is Release)" FORCE)
endif()

# Set parallel build level if not already defined
if(NOT DEFINED ENV{CMAKE_BUILD_PARALLEL_LEVEL})
    # Detect the number of physical cores instead of logical threads
    execute_process(
        COMMAND bash -c "lscpu | grep '^Core(s) per socket:' | awk '{print $4}'"
        OUTPUT_VARIABLE CORES_PER_SOCKET OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    execute_process(
        COMMAND bash -c "lscpu | grep '^Socket(s):' | awk '{print $2}'"
        OUTPUT_VARIABLE SOCKETS OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    math(EXPR PHYSICAL_CORES "${CORES_PER_SOCKET} * ${SOCKETS}")
    # Set CMAKE_BUILD_PARALLEL_LEVEL
    set(ENV{CMAKE_BUILD_PARALLEL_LEVEL} ${PHYSICAL_CORES})
    message(STATUS "CMAKE_BUILD_PARALLEL_LEVEL set to ${PHYSICAL_CORES} (physical cores)")
else()
    message(STATUS "CMAKE_BUILD_PARALLEL_LEVEL is already set to $ENV{CMAKE_BUILD_PARALLEL_LEVEL}")
endif()

# Set the install directory for GSL
if (NOT GSL_ROOT_DIR)
    set(GSL_ROOT_DIR ${CMAKE_BINARY_DIR}/gsl-install)
endif()

# Find GSL library
message(STATUS "Searching for GSL version>=2.8 at ${GSL_ROOT_DIR}")

find_package(GSL 2.8 MODULE)

if (NOT GSL_FOUND)
    message(STATUS "GSL version>=2.8 not found (${GSL_ROOT_DIR}). Downloading and building GSL")
    include(FetchContent)

    # Set install location for deal.II
    # Fetch deal.II
    FetchContent_Declare(
        gsl
        URL https://ftp.gnu.org/gnu/gsl/gsl-latest.tar.gz
        SOURCE_DIR ${CMAKE_BINARY_DIR}/gsl-src
        INSTALL_DIR ${GSL_ROOT_DIR}
    )

    # Ensure deal.II is built before calling find_package
    FetchContent_GetProperties(gsl)
    if(NOT gsl_POPULATED)
        FetchContent_Populate(gsl)
        execute_process(
            COMMAND ./configure --prefix=${GSL_ROOT_DIR} WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/gsl-src
        )
        # Build deal.II
        execute_process(
            COMMAND make -j$ENV{CMAKE_BUILD_PARALLEL_LEVEL} WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/gsl-src RESULT_VARIABLE MAKE_RESULT
        )

        if(NOT MAKE_RESULT EQUAL 0)
            message(FATAL_ERROR "Failed to build GSL with make")
        endif()

        # Install deal.II
        execute_process(
            COMMAND make install WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/gsl-src
        )
    endif()

    message(STATUS "GSL installed to ${GSL_ROOT_DIR}")
    find_package(GSL 2.8 REQUIRED MODULE)
endif()

# Find TBB library
find_package(TBB REQUIRED)

# Include directories
include_directories(${PROJECT_SOURCE_DIR}/include)

# Add the source files
file(GLOB SOURCES ${PROJECT_SOURCE_DIR}/src/*.cpp)

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

if(ENABLE_JAVA_WRAPPER)
    # Find Java only if JNI is enabled
    find_package(Java REQUIRED)
    find_package(JNI REQUIRED)

    include(UseJava)  # CMake module to compile Java
    
    # Add the source files
    file(GLOB_RECURSE JAVA_SOURCES ${PROJECT_SOURCE_DIR}/src/jni/*.java)  # Add Java interface

    # Compile Java files
    add_jar(GetelecJNI ${JAVA_SOURCES} OUTPUT_NAME GetelecJNI)


    # Include JNI headers
    include_directories(${JNI_INCLUDE_DIRS})

    # Add JNI shared library target
    list(APPEND SOURCES ${PROJECT_SOURCE_DIR}/src/jni/*.cpp)


    # Link JNI libraries if needed
    target_link_libraries(getelec_lib ${JNI_LIBRARIES})

    add_custom_command(
        TARGET GetelecJNI POST_BUILD
        COMMAND ${Java_JAVAC_EXECUTABLE} -cp ${CMAKE_BINARY_DIR}/GetelecJNI.jar -d ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/src/jni/Main.java
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        COMMENT "Compiling Main.java into build directory"
    )

    message(STATUS "JNI support enabled. Compiling Java wrapper.")
endif()

# Add debug flags if building in Debug mode
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_definitions(getelec_lib PRIVATE DEBUG_MODE)
    target_compile_options(getelec_lib PRIVATE -g -Wall -Wextra -Wno-reorder -Wno-unused-value -O0 -DDEBUG)
    target_compile_definitions(getelec PRIVATE DEBUG_MODE)
    target_compile_options(getelec PRIVATE -g -O0 -DDEBUG)
endif()


# Add release-specific flags if building in Release mode
if (CMAKE_BUILD_TYPE STREQUAL "Release")
    target_compile_options(getelec_lib PRIVATE -O3 -DRELEASE)
    target_compile_options(getelec PRIVATE -O3 -DRELEASE)
endif()


if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT OR NOT DEFINED CMAKE_INSTALL_PREFIX)
    set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/install" CACHE PATH "Install path prefix" FORCE)
    message(STATUS "Setting default CMAKE_INSTALL_PREFIX path to ${CMAKE_INSTALL_PREFIX}")
endif()


# Install the getelec library locally
install(TARGETS getelec_lib EXPORT GetelecTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
    ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
    RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
    INCLUDES DESTINATION ${CMAKE_INSTALL_PREFIX}/include
)

# Install the getelec executable
install(TARGETS getelec
    RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
)

# # Install include headers
install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/ DESTINATION include)


# Export the GetelecConfig.cmake for find_package()
install(EXPORT GetelecTargets
    FILE GetelecConfig.cmake
    NAMESPACE getelec::
    DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/cmake/getelec
)

# Allow this project to be used with `find_package(getelec REQUIRED)`
export(EXPORT GetelecTargets FILE ${CMAKE_BINARY_DIR}/GetelecConfig.cmake NAMESPACE getelec::)

# Append find_dependency() calls directly to the installed GetelecConfig.cmake file
install(CODE "
    file(APPEND \"\${CMAKE_INSTALL_PREFIX}/lib/cmake/getelec/GetelecConfig.cmake\" \"\\ninclude(CMakeFindDependencyMacro)\")
    file(APPEND \"\${CMAKE_INSTALL_PREFIX}/lib/cmake/getelec/GetelecConfig.cmake\" \"\\nfind_dependency(TBB REQUIRED)\")
    file(APPEND \"\${CMAKE_INSTALL_PREFIX}/lib/cmake/getelec/GetelecConfig.cmake\" \"\\nfind_dependency(GSL REQUIRED)\")
")

option(BUILD_TESTS "Build tests" OFF)

if (BUILD_TESTS)

    #enable the tests
    enable_testing()

    # Include the tests directory
    add_subdirectory(tests)
endif()

