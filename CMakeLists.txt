cmake_minimum_required(VERSION 3.9.1)
project(k_means)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add OpenMP support
find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

set(GCC_COVERAGE_COMPILE_FLAGS "-pedantic")

# Add the include directories
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

# Source files (ensure all are added)
file(GLOB SOURCES CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)

# Create the executable
add_executable(k_means ${SOURCES})

# Link OpenMP and MPI if necessary
find_package(MPI REQUIRED)
target_link_libraries(k_means PUBLIC MPI::MPI_CXX)

# Link OpenMP explicitly (if not automatically included)
if(OpenMP_CXX_FOUND)
    target_link_libraries(k_means PUBLIC OpenMP::OpenMP_CXX)
endif()

target_link_libraries(k_means PUBLIC stdc++)

# DEBUG Options
#target_compile_options(k_means PRIVATE -fsanitize=address)
#target_link_options(k_means PRIVATE -fsanitize=address)

