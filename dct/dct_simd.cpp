/************************************************************************************
 **  Highway SIMD DCT Implementation                                              **
 **  Portable SIMD that works on x86, ARM, RISC-V, etc.                          **
 ************************************************************************************/

#include "dct/dct_simd.hpp"
#include <cstring>

// Highway for portable SIMD
#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "dct/dct_simd.cpp"
#include <hwy/foreach_target.h>
#include <hwy/highway.h>

HWY_BEFORE_NAMESPACE();
namespace dct_simd {
namespace HWY_NAMESPACE {

namespace hn = hwy::HWY_NAMESPACE;

// Hamming window coefficients (precomputed)
// cos((x-3.5)/4.5 * pi/2) * cos((y-3.5)/4.5 * pi/2) / (norm/64)
static const float window_table[8][8] = {
  {0.2328f, 0.4375f, 0.5894f, 0.6702f, 0.6702f, 0.5894f, 0.4375f, 0.2328f},
  {0.4375f, 0.8222f, 1.1077f, 1.2596f, 1.2596f, 1.1077f, 0.8222f, 0.4375f},
  {0.5894f, 1.1077f, 1.4924f, 1.6971f, 1.6971f, 1.4924f, 1.1077f, 0.5894f},
  {0.6702f, 1.2596f, 1.6971f, 1.9298f, 1.9298f, 1.6971f, 1.2596f, 0.6702f},
  {0.6702f, 1.2596f, 1.6971f, 1.9298f, 1.9298f, 1.6971f, 1.2596f, 0.6702f},
  {0.5894f, 1.1077f, 1.4924f, 1.6971f, 1.6971f, 1.4924f, 1.1077f, 0.5894f},
  {0.4375f, 0.8222f, 1.1077f, 1.2596f, 1.2596f, 1.1077f, 0.8222f, 0.4375f},
  {0.2328f, 0.4375f, 0.5894f, 0.6702f, 0.6702f, 0.5894f, 0.4375f, 0.2328f}
};

// Highway version of 1D DCT for 8 values
// Uses Highway for initial butterfly, scalar for DCT butterfly (which is already scalar-heavy)
void DCT_1D_Row_HWY(const float* __restrict__ row, float* __restrict__ out) {
  // Use scalable tag - Highway will use appropriate vector size
  const hn::ScalableTag<float> d;
  
  // For 8 values, we'll process them (may be less than 8 if vector is smaller)
  // Load up to N floats
  auto vec = hn::LoadU(d, row);
  
  // Reverse the vector for butterfly
  auto reversed = hn::Reverse(d, vec);
  
  // Butterfly: sum and difference
  auto sum_vec = hn::Add(vec, reversed);
  auto diff_vec = hn::Sub(vec, reversed);
  
  // Store to arrays for scalar DCT butterfly
  alignas(64) float tmp_sum[8], tmp_diff[8];
  hn::StoreU(sum_vec, d, tmp_sum);
  hn::StoreU(diff_vec, d, tmp_diff);
  
  // Extract values (sum_vec has symmetry: [tmp0 tmp1 tmp2 tmp3 tmp3 tmp2 tmp1 tmp0])
  float tmp0 = tmp_sum[0], tmp1 = tmp_sum[1], tmp2 = tmp_sum[2], tmp3 = tmp_sum[3];
  float tmp7 = tmp_diff[0], tmp6 = tmp_diff[1], tmp5 = tmp_diff[2], tmp4 = tmp_diff[3];
  
  // Even part (DCT butterfly - inherently scalar due to dependencies)
  float tmp10 = tmp0 + tmp3;
  float tmp13 = tmp0 - tmp3;
  float tmp11 = tmp1 + tmp2;
  float tmp12 = tmp1 - tmp2;
  
  float out0 = tmp10 + tmp11;
  float out4 = tmp10 - tmp11;
  
  float z1 = (tmp12 + tmp13) * 0.707106781f;  // 1/sqrt(2)
  float out2 = tmp13 + z1;
  float out6 = tmp13 - z1;
  
  // Odd part
  tmp10 = tmp4 + tmp5;
  tmp11 = tmp5 + tmp6;
  tmp12 = tmp6 + tmp7;
  
  float z5 = (tmp10 - tmp12) * 0.382683433f;  // sqrt(2) - 1
  float z2 = 0.541196100f * tmp10 + z5;       // cos(3*pi/8)
  float z4 = 1.306562965f * tmp12 + z5;       // cos(pi/8)
  float z3 = tmp11 * 0.707106781f;
  
  float z11 = tmp7 + z3;
  float z13 = tmp7 - z3;
  
  float out5 = z13 + z2;
  float out3 = z13 - z2;
  float out1 = z11 + z4;
  float out7 = z11 - z4;
  
  // Store results
  out[0] = out0; out[1] = out1; out[2] = out2; out[3] = out3;
  out[4] = out4; out[5] = out5; out[6] = out6; out[7] = out7;
}

// 2D DCT: Apply 1D DCT to rows, transpose, apply 1D DCT to rows (columns of original), apply window
void DCT_2D_Block_HWY(const float* __restrict__ block, float* __restrict__ output) {
  // Temporary storage for row DCT and transpose
  alignas(64) float row_dct[64];
  alignas(64) float transposed[64];
  
  // Step 1: Apply 1D DCT to each row
  for (int y = 0; y < 8; y++) {
    DCT_1D_Row_HWY(&block[y * 8], &row_dct[y * 8]);
  }
  
  // Step 2: Transpose (scalar - simple and cache-friendly)
  for (int y = 0; y < 8; y++) {
    for (int x = 0; x < 8; x++) {
      transposed[x * 8 + y] = row_dct[y * 8 + x];
    }
  }
  
  // Step 3: Apply 1D DCT to each row (which are columns of original)
  alignas(64) float col_dct[64];
  for (int y = 0; y < 8; y++) {
    DCT_1D_Row_HWY(&transposed[y * 8], &col_dct[y * 8]);
  }
  
  // Step 4: Apply Hamming window and transpose back
  const hn::ScalableTag<float> d;
  for (int y = 0; y < 8; y++) {
    // Load row of DCT coefficients
    auto dct_row = hn::LoadU(d, &col_dct[y * 8]);
    
    // Load window coefficients
    auto window = hn::LoadU(d, window_table[y]);
    
    // Multiply
    auto result = hn::Mul(dct_row, window);
    
    // Transpose while storing (store column-wise)
    alignas(64) float temp[8] = {0};  // Initialize temp array
    hn::StoreU(result, d, temp);
    for (int x = 0; x < 8; x++) {
      output[x * 8 + y] = temp[x];
    }
  }
}

float GetWindowCoeff_HWY(int x, int y) {
  return window_table[y][x];
}

} // namespace HWY_NAMESPACE
} // namespace dct_simd
HWY_AFTER_NAMESPACE();

#if HWY_ONCE
namespace dct_simd {

// Export functions for dynamic dispatch
HWY_EXPORT(DCT_1D_Row_HWY);
HWY_EXPORT(DCT_2D_Block_HWY);
HWY_EXPORT(GetWindowCoeff_HWY);

// Helper function to get all compiled Highway targets from this multi-target module
// This works because this file is compiled with foreach_target.h
std::vector<int64_t> GetCompiledTargets() {
  std::vector<int64_t> targets;
  for (const int64_t target : hwy::SupportedAndGeneratedTargets()) {
    targets.push_back(target);
  }
  return targets;
}

// Wrapper functions
void DCT_1D_Row(const float* row, float* out) {
  HWY_DYNAMIC_DISPATCH(DCT_1D_Row_HWY)(row, out);
}

void DCT_2D_Block(const float* block, float* output) {
  HWY_DYNAMIC_DISPATCH(DCT_2D_Block_HWY)(block, output);
}

float GetWindowCoeff(int x, int y) {
  return HWY_DYNAMIC_DISPATCH(GetWindowCoeff_HWY)(x, y);
}

} // namespace dct_simd
#endif // HWY_ONCE
