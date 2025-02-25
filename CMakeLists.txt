# Set the minimum version of CMake required
cmake_minimum_required(VERSION 3.10)

# Define the project name and version
project(GETELEC VERSION 3.0 LANGUAGES CXX)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

option(ENABLE_JAVA_WRAPPER "Enable Java JNI Wrapper" OFF)


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
    target_compile_options(getelec_lib PRIVATE -g -O0 -DDEBUG)
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

# Export the GetelecConfig.cmake for find_package()
install(EXPORT GetelecTargets
    FILE GetelecConfig.cmake
    NAMESPACE getelec::
    DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/cmake/getelec
)


# # Install include headers
# install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/ DESTINATION include)

# Export the GetelecConfig.cmake for find_package()
install(EXPORT GetelecTargets
    FILE GetelecConfig.cmake
    NAMESPACE getelec::
    DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/cmake/getelec
)

# Allow this project to be used with `find_package(getelec REQUIRED)`
export(EXPORT GetelecTargets FILE ${CMAKE_BINARY_DIR}/GetelecConfig.cmake NAMESPACE getelec::)

option(BUILD_TESTS "Build tests" OFF)

if (BUILD_TESTS)

    #enable the tests
    enable_testing()

    # Include the tests directory
    add_subdirectory(tests)
endif()

