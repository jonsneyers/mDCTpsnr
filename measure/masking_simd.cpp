/************************************************************************************
 **  Highway SIMD implementations for masking operations                         **
 **  Portable SIMD that works on x86 (SSE2/AVX2/AVX-512), ARM (NEON), etc.       **
 ************************************************************************************/

#include "measure/masking_simd.hpp"

// Highway for portable SIMD
#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "measure/masking_simd.cpp"
#include <hwy/foreach_target.h>
#include <hwy/highway.h>

HWY_BEFORE_NAMESPACE();
namespace masking_simd {
namespace HWY_NAMESPACE {

namespace hn = hwy::HWY_NAMESPACE;

// Highway version of slope=1.0 loop (fabsf + multiply)
__attribute__((hot))
void ApplyMaskingSlope1_HWY(const float* __restrict__ input_ptr, 
                             float* __restrict__ mapped_ptr, 
                             float vis, size_t width) {
  const hn::ScalableTag<float> d;
  const size_t N = hn::Lanes(d);
  auto vvis = hn::Set(d, vis);
  
  // Buffers are 64-byte aligned, use aligned loads/stores for better performance
  #pragma GCC ivdep
  
  for (size_t i = 0; i < width; i += N) {
    auto vin = hn::Load(d, &input_ptr[i]);
    auto vabs = hn::Abs(vin);
    auto vout = hn::Mul(vabs, vvis);
    hn::Store(vout, d, &mapped_ptr[i]);
  }
}

// Highway version of slope=0.5 loop (sqrt of scaled absolute)
__attribute__((hot))
void ApplyMaskingSlope05_HWY(const float* __restrict__ input_ptr, 
                              float* __restrict__ mapped_ptr,
                              float vis, size_t width) {
  const hn::ScalableTag<float> d;
  const size_t N = hn::Lanes(d);
  auto vvis = hn::Set(d, vis);
  
  // Buffers are 64-byte aligned, use aligned loads/stores
  #pragma GCC ivdep
  
  for (size_t i = 0; i < width; i += N) {
    auto vin = hn::Load(d, &input_ptr[i]);
    auto vabs = hn::Abs(vin);
    auto vmul = hn::Mul(vabs, vvis);
    auto vout = hn::Sqrt(vmul);
    hn::Store(vout, d, &mapped_ptr[i]);
  }
}

// Highway version of column sum convolution (5-tap symmetric filter)
__attribute__((hot))
void ColumnSumConvolution_HWY(const float* __restrict__ mapped_ptr, 
                               float* __restrict__ added_ptr,
                               float m0, float m1, float m2, size_t width) {
  const hn::ScalableTag<float> d;
  const size_t N = hn::Lanes(d);
  auto vm0 = hn::Set(d, m0);
  auto vm1 = hn::Set(d, m1);
  auto vm2 = hn::Set(d, m2);
  
  // Buffers are 64-byte aligned, use aligned loads/stores
  // Note: We load from offset positions which may not be aligned
  #pragma GCC ivdep
  
  for (size_t i = 0; i < width; i += N) {
    // Exploit symmetry: m0*(v0+v4) + m1*(v1+v3) + m2*v2
    auto v0 = hn::LoadU(d, &mapped_ptr[i]);      // LoadU needed - offset 0
    auto v1 = hn::LoadU(d, &mapped_ptr[i + 1]);  // LoadU needed - offset 1
    auto v2 = hn::LoadU(d, &mapped_ptr[i + 2]);  // LoadU needed - offset 2
    auto v3 = hn::LoadU(d, &mapped_ptr[i + 3]);  // LoadU needed - offset 3
    auto v4 = hn::LoadU(d, &mapped_ptr[i + 4]);  // LoadU needed - offset 4
    
    auto sum04 = hn::Add(v0, v4);
    auto sum13 = hn::Add(v1, v3);
    
    auto term0 = hn::Mul(vm0, sum04);
    auto term1 = hn::Mul(vm1, sum13);
    auto term2 = hn::Mul(vm2, v2);
    
    auto result = hn::Add(hn::Add(term0, term1), term2);
    hn::Store(result, d, &added_ptr[i]);  // Output is aligned
  }
}

// Highway version of ComputeMask main loop
__attribute__((hot))
void ComputeMaskLoop_HWY(float* __restrict__ * __restrict__ added_ptrs, 
                         const float* __restrict__ input_center,
                         float* __restrict__ mask_out, 
                         float* __restrict__ output_out,
                         float m0, float m1, float m2,
                         float base, float invMaskSizeSq, float visibility,
                         size_t width, bool no_base_visibility) {
  const hn::ScalableTag<float> d;
  const size_t N = hn::Lanes(d);
  auto vm0 = hn::Set(d, m0);
  auto vm1 = hn::Set(d, m1);
  auto vm2 = hn::Set(d, m2);
  auto vbase = hn::Set(d, base);
  auto vinvMaskSizeSq = hn::Set(d, invMaskSizeSq);
  auto vvis = hn::Set(d, visibility);
  auto vone = hn::Set(d, 1.0f);
  
  // Buffers are 64-byte aligned, use aligned loads/stores
  #pragma GCC ivdep
  
  for (size_t i = 0; i < width; i += N) {
    // Load 5 rows of added data (exploiting symmetry)
    auto vdata0 = hn::Load(d, &added_ptrs[0][i]);
    auto vdata1 = hn::Load(d, &added_ptrs[1][i]);
    auto vdata2 = hn::Load(d, &added_ptrs[2][i]);
    auto vdata3 = hn::Load(d, &added_ptrs[3][i]);
    auto vdata4 = hn::Load(d, &added_ptrs[4][i]);
    
    // Compute weighted sum exploiting symmetry: m0*(d0+d4) + m1*(d1+d3) + m2*d2
    auto sum04 = hn::Add(vdata0, vdata4);
    auto sum13 = hn::Add(vdata1, vdata3);
    
    auto vsum = hn::Mul(vm0, sum04);
    vsum = hn::MulAdd(vm1, sum13, vsum);
    vsum = hn::MulAdd(vm2, vdata2, vsum);
    
    vsum = hn::Mul(vsum, vinvMaskSizeSq);
    
    // Compute mask result
    auto vmask_result = no_base_visibility ? 
      hn::Div(vone, vsum) :
      hn::Div(vbase, hn::Add(vbase, vsum));
    
    hn::Store(vmask_result, d, &mask_out[i]);
    
    // Compute output
    auto vinput = hn::Load(d, &input_center[i]);
    auto voutput = hn::Mul(vinput, vvis);
    hn::Store(voutput, d, &output_out[i]);
  }
}

}  // namespace HWY_NAMESPACE
}  // namespace masking_simd
HWY_AFTER_NAMESPACE();

#if HWY_ONCE
namespace masking_simd {

// Export functions for dynamic dispatch
HWY_EXPORT(ApplyMaskingSlope1_HWY);
HWY_EXPORT(ApplyMaskingSlope05_HWY);
HWY_EXPORT(ColumnSumConvolution_HWY);
HWY_EXPORT(ComputeMaskLoop_HWY);

// Wrapper functions
__attribute__((hot, flatten))
void ApplyMaskingSlope1(const float* __restrict__ input_ptr, 
                        float* __restrict__ mapped_ptr,
                        float vis, size_t width) {
  return HWY_DYNAMIC_DISPATCH(ApplyMaskingSlope1_HWY)(input_ptr, mapped_ptr, vis, width);
}

__attribute__((hot, flatten))
void ApplyMaskingSlope05(const float* __restrict__ input_ptr, 
                         float* __restrict__ mapped_ptr,
                         float vis, size_t width) {
  return HWY_DYNAMIC_DISPATCH(ApplyMaskingSlope05_HWY)(input_ptr, mapped_ptr, vis, width);
}

__attribute__((hot, flatten))
void ColumnSumConvolution(const float* __restrict__ mapped_ptr, 
                          float* __restrict__ added_ptr,
                          float m0, float m1, float m2, size_t width) {
  return HWY_DYNAMIC_DISPATCH(ColumnSumConvolution_HWY)(mapped_ptr, added_ptr, m0, m1, m2, width);
}

__attribute__((hot, flatten))
void ComputeMaskLoop(float* __restrict__ * __restrict__ added_ptrs, 
                     const float* __restrict__ input_center,
                     float* __restrict__ mask_out, 
                     float* __restrict__ output_out,
                     float m0, float m1, float m2,
                     float base, float invMaskSizeSq, float visibility,
                     size_t width, bool no_base_visibility) {
  return HWY_DYNAMIC_DISPATCH(ComputeMaskLoop_HWY)(added_ptrs, input_center, mask_out, output_out,
                                                     m0, m1, m2, base, invMaskSizeSq, visibility,
                                                     width, no_base_visibility);
}

} // namespace masking_simd
#endif // HWY_ONCE
