/************************************************************************************
 **  Highway SIMD functions for masking operations                               **
 **  Portable SIMD that works on x86 (SSE2/AVX2/AVX-512), ARM (NEON), etc.       **
 ************************************************************************************/

#ifndef MASKING_SIMD_HPP
#define MASKING_SIMD_HPP

#include <stddef.h>

namespace masking_simd {

// Highway version of slope=1.0 loop (fabsf + multiply)
void ApplyMaskingSlope1(const float* __restrict__ input_ptr, 
                        float* __restrict__ mapped_ptr,
                        float vis, size_t width);

// Highway version of slope=0.5 loop (sqrt of scaled absolute)
void ApplyMaskingSlope05(const float* __restrict__ input_ptr, 
                         float* __restrict__ mapped_ptr,
                         float vis, size_t width);

// Highway version of column sum convolution (5-tap symmetric filter)
void ColumnSumConvolution(const float* __restrict__ mapped_ptr, 
                          float* __restrict__ added_ptr,
                          float m0, float m1, float m2, size_t width);

// Highway version of ComputeMask main loop
void ComputeMaskLoop(float* __restrict__ * __restrict__ added_ptrs, 
                     const float* __restrict__ input_center, 
                     float* __restrict__ mask_out, 
                     float* __restrict__ output_out,
                     float m0, float m1, float m2,
                     float base, float invMaskSizeSq, float visibility,
                     size_t width, bool no_base_visibility);

}  // namespace masking_simd

#endif  // MASKING_SIMD_HPP

