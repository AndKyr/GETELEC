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

# Install the getelec library
install(TARGETS getelec_lib EXPORT GetelecTargets
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    RUNTIME DESTINATION bin
    INCLUDES DESTINATION include
)

# Install include headers
install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/ DESTINATION include)

# Export the getelec target for other projects to use
install(EXPORT GetelecTargets
    FILE GetelecConfig.cmake
    NAMESPACE getelec::
    DESTINATION lib/cmake/getelec
)

# Allow this project to be used with `find_package(getelec REQUIRED)`
export(EXPORT GetelecTargets FILE ${CMAKE_BINARY_DIR}/GetelecConfig.cmake NAMESPACE getelec::)



#enable the tests
enable_testing()

# Include the tests directory
add_subdirectory(tests)

