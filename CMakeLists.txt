cmake_minimum_required(VERSION 3.15)
project(dctpsnr VERSION 1.0 LANGUAGES CXX)

# Use Clang by default
if(NOT CMAKE_CXX_COMPILER)
    find_program(CLANG_CXX clang++-18 clang++)
    if(CLANG_CXX)
        set(CMAKE_CXX_COMPILER ${CLANG_CXX})
    endif()
endif()

# C++11 required
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wunused -Wpointer-arith -pedantic")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wcast-qual -Wwrite-strings -Wno-long-long -Winline")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wredundant-decls -Wnon-virtual-dtor -Woverloaded-virtual")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wsign-promo -Wshadow -Wuninitialized")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wdouble-promotion -Wfloat-conversion")

# Optimization flags for Release
# Note: We don't use -march=native globally because Highway needs to compile
# for multiple targets. Instead, we'll add it to non-Highway files.
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -mtune=native -ffast-math")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -funroll-loops -fstrict-aliasing")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -ftree-vectorize -fomit-frame-pointer")
# Additional math optimizations
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -fno-math-errno -fno-trapping-math")
# Better inlining
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -finline-limit=1000")
# Note: LTO disabled - causes performance regression with Highway's multi-target compilation

# Debug flags
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g3 -fno-inline -DDEBUG_LEVEL=2")

# Find Highway library
find_package(PkgConfig REQUIRED)
pkg_check_modules(HWY REQUIRED libhwy)

# Include directories
include_directories(${CMAKE_SOURCE_DIR})
include_directories(${HWY_INCLUDE_DIRS})

# Highway-specific flags
add_definitions(-DHWY_SHARED_DEFINE)

# Source files organized by directory
set(CMD_SOURCES
    cmd/main.cpp
)

set(GLOBAL_SOURCES
    global/exceptions.cpp
    global/thread.cpp
)

set(IO_SOURCES
    io/bytestream.cpp
    io/filestream.cpp
)

set(IMG_SOURCES
    img/image.cpp
    img/imgwriter.cpp
)

set(CTRAFO_SOURCES
    ctrafo/colortransformer.cpp
)

set(DCT_SOURCES
    dct/line.cpp
    dct/component.cpp
    dct/dct_simd.cpp
)

set(MEASURE_SOURCES
    measure/masking.cpp
    measure/masking_simd.cpp
    measure/pooling.cpp
    measure/ediff.cpp
)

# Create executable
add_executable(dctpsnr
    ${CMD_SOURCES}
    ${GLOBAL_SOURCES}
    ${IO_SOURCES}
    ${IMG_SOURCES}
    ${CTRAFO_SOURCES}
    ${DCT_SOURCES}
    ${MEASURE_SOURCES}
)

# For Highway to compile multiple targets, we need to let it use its own target selection
# Don't use -march=native on Highway SIMD files, but use it for other files
set_source_files_properties(
    ${CMD_SOURCES}
    ${GLOBAL_SOURCES}
    ${IO_SOURCES}
    ${IMG_SOURCES}
    ${CTRAFO_SOURCES}
    ${DCT_SOURCES}
    measure/masking.cpp
    measure/pooling.cpp
    measure/ediff.cpp
    PROPERTIES COMPILE_FLAGS "-march=native"
)

# For Highway SIMD files, compile for all attainable targets (SSE4, AVX2, AVX-512, etc.)
# This enables dynamic dispatch to select the best available target at runtime
# Use HWY_IS_TEST and ensure HWY_SKIP_NON_BEST_BASELINE is not defined
set_source_files_properties(
    measure/masking_simd.cpp
    measure/pooling.cpp
    dct/dct_simd.cpp
    dct/component.cpp
    PROPERTIES COMPILE_DEFINITIONS "HWY_IS_TEST;HWY_SKIP_NON_BEST_BASELINE=0" 
               COMPILE_FLAGS "-mavx512f -mavx512vl -mavx512dq -mavx512bw -mavx2 -mfma -mavx -msse4.2 -mssse3"
)

# Link libraries
target_link_libraries(dctpsnr
    ${HWY_LIBRARIES}
    m
)

# Installation
install(TARGETS dctpsnr DESTINATION bin)

# Print configuration
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "C++ compiler: ${CMAKE_CXX_COMPILER}")
message(STATUS "C++ flags: ${CMAKE_CXX_FLAGS}")
message(STATUS "C++ flags (Release): ${CMAKE_CXX_FLAGS_RELEASE}")
message(STATUS "Highway found: ${HWY_FOUND}")
message(STATUS "Highway include dirs: ${HWY_INCLUDE_DIRS}")
message(STATUS "Highway libraries: ${HWY_LIBRARIES}")

