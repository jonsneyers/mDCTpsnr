# mDCTPSNR - Modified DCT-based Perceptual Similarity Metric

A high-performance, portable SIMD implementation of the modified DCT-based perceptual image quality metric for measuring just-noticeable differences (JND) in compressed images.

## Overview

mDCTPSNR is a full-reference image quality metric that operates in the DCT (Discrete Cosine Transform) domain to measure perceptual differences between a reference image and a distorted version. Unlike traditional PSNR, it incorporates:

- **Perceptual masking**: Frequency-dependent visibility thresholds based on human visual system models
- **Local error pooling**: Spatial aggregation that accounts for the non-uniform distribution of errors
- **DCT-domain operation**: Efficient computation aligned with common image compression formats

The metric outputs a JND (Just Noticeable Difference) score where:
- **Lower values** = higher quality (less perceptible difference)
- **~1.0** = threshold of visibility  
- **Higher values** = more visible distortion

## Key Features

- **Portable SIMD**: Uses [Google Highway](https://github.com/google/highway) for cross-platform vectorization
  - Automatically selects best available: AVX-512, AVX2, SSE4, ARM NEON, etc.
  - Runtime CPU dispatch - single binary works optimally on any processor
  - No manual intrinsics - Highway abstracts platform differences
  
- **High Performance**: ~7x faster than original implementation
  - Optimized hot paths with SIMD vectorization
  - Pre-computed lookup tables for non-linear transformations
  - Efficient memory layout and cache utilization
  
- **Modern Codebase**: 
  - CMake build system with Clang 18
  - C++11 standard library types
  - Clean architecture without legacy code paths

## Building

### Requirements

- **Compiler**: Clang 18+ (recommended) or GCC 11+
- **CMake**: 3.15 or later
- **Highway**: libhwy 1.2.0 or later (`sudo apt install libhwy-dev` or `sudo dnf install highway-devel`)

### Compilation

```bash
git clone https://github.com/thorfdbg/mDCTpsnr.git
cd mDCTpsnr

# Build
mkdir build && cd build
# CC=clang-18 CXX=clang++-18
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build . -j$(nproc)

# Binary location: build/dctpsnr
```

## Usage

```bash
# Convert images to PPM format if needed (required input format)
convert reference.png reference.ppm
convert distorted.png distorted.ppm

# Compute JND score
./dctpsnr -jnd reference.ppm distorted.ppm

# Output: single floating-point JND value
# Example: 2.54 (means ~2.5 JND units of visible distortion)
```

### Multi-threading

```bash
# Use 4 threads for faster processing on large images
./dctpsnr -jnd -P 4 reference.ppm distorted.ppm
```

### Version Information

```bash
# Show compiled SIMD targets and CPU capabilities
./dctpsnr --version
```

## Algorithm

### Processing Pipeline

1. **Image Loading**: Read PPM files, store raw integer pixel values
2. **Color Transform**: Convert sRGB → Linear RGB via pre-computed lookup table, then convert to Linear YCbCr
3. **DCT Transform**: Apply 8x8 DCT to both reference and distorted images
4. **Frequency Masking**: Apply perceptual visibility thresholds per DCT coefficient
5. **Error Computation**: Calculate masked error between DCT coefficients
6. **Spatial Pooling**: Aggregate errors spatially using Minkowski summation
7. **JND Output**: Final score representing perceptual difference

### Key Parameters

Configured in `options.hpp`:

- **Color Space**: `LINEAR` (Linear RGB, default) or `OPPONENT_COLOR`
- **Masking Model**: Ahumada model with exponent 3.5
- **Base Visibility**: 0.08 (threshold at which artifacts become visible)
- **Detection Threshold**: 1.0 JND


## Architecture

```
dctpsnr-highway/
├── cmd/              Main entry point and command-line parsing
├── global/           Core utilities, threading, exceptions
├── io/               File I/O (ByteStream abstraction)
├── img/              Image loading and buffer management
├── ctrafo/           Color space transformations (sRGB ↔ Linear)
├── dct/              DCT operations, line buffers, SIMD transforms
├── measure/          Core metric computation
│   ├── masking.cpp          Perceptual masking functions
│   ├── masking_simd.cpp     Highway SIMD implementations
│   ├── pooling.cpp          Error aggregation and JND computation
│   └── ediff.cpp            Error difference calculation
├── test/             Test images, scripts, and baselines
└── docs-optimization/  Detailed optimization documentation
```

### SIMD Optimization

Critical loops are vectorized using Google Highway:

- **Masking operations**: 4-16x parallelism depending on CPU vector width
- **sRGB linearization**: O(1) lookup table replaces expensive `powf()` calls
- **DCT transforms**: Portable 8x8 matrix operations
- **Error pooling**: Vectorized aggregation

Highway provides a single source implementation that compiles to:
- **AVX-512**: 16 floats/vector (modern Intel/AMD)
- **AVX2**: 8 floats/vector (Intel Haswell+, AMD Zen+)
- **SSE4**: 4 floats/vector (legacy x86)
- **NEON**: 4 floats/vector (ARM)


## License

Original implementation copyright © 2009 University of Stuttgart, Thomas Richter.
Provided 'as-is' with permission for use and modification (see source file headers).

Highway port and optimizations © 2026 Cloudinary, Jon Sneyers.
Provided 'as-is' with permission for use and modification (see source file headers).


