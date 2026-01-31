/************************************************************************************
 **  Highway SIMD DCT operations - Portable SIMD for DCT computation             **
 **  Works on x86 (SSE2/AVX2/AVX-512), ARM (NEON), RISC-V, etc.                  **
 ************************************************************************************/

#ifndef DCT_DCT_SIMD_HPP
#define DCT_DCT_SIMD_HPP

#include <cstddef>
#include <vector>
#include <cstdint>

namespace dct_simd {

// Get Hamming window coefficient for given position
float GetWindowCoeff(int x, int y);

// 1D DCT for a row of 8 floats
// Input: 8 consecutive floats
// Output: 8 DCT coefficients (writes to out array)
void DCT_1D_Row(const float* row, float* out);

// Full 2D DCT for an 8x8 block
// Input: 8x8 block (row-major, 64 floats)
// Output: 8x8 DCT coefficients with Hamming window applied
void DCT_2D_Block(const float* block, float* output);

// Get all compiled Highway targets from this multi-target module
std::vector<int64_t> GetCompiledTargets();

} // namespace dct_simd

#endif // DCT_DCT_SIMD_HPP
