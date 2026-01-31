/************************************************************************************
 **  Copyright (C) 2009 University of Stuttgart, Thomas Richter                    **
 **                                                                                **
 **  This software is provided 'as-is', without any express or implied             **
 **  warranty.  In no event will the authors be held liable for any damages        **
 **  arising from the use of this software.                                        **
 **                                                                                **
 **  Permission is granted to anyone to use this software for any purpose,         **
 **  including commercial applications, and to alter it and redistribute it        **
 **  freely, subject to the following restrictions:                                **
 **                                                                                **
 **  1. The origin of this software must not be misrepresented; you must not       **
 **     claim that you wrote the original software. If you use this software       **
 **     in a product, an acknowledgment in the product documentation would be      **
 **     appreciated but is not required.                                           **
 **  2. Altered source versions must be plainly marked as such, and must not be    **
 **     misrepresented as being the original software.                             **
 **  3. This notice may not be removed or altered from any source distribution.    **
 **                                                                                **
 **     Thomas Richter                                                             **
 **     thomas.richter@iis.fraunhofer.de                                           **
 **                                                                                **
 ************************************************************************************/

/*
** $Id: component.cpp,v 1.9 2025/07/08 14:35:22 thor Exp $
**
*/

/// Includes
#include <cstdint>
#include <cassert>
#include <cmath>
#include "dct/line.hpp"
#include "dct/component.hpp"

// Highway for portable SIMD
#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "dct/component.cpp"
#include <hwy/foreach_target.h>
#include <hwy/highway.h>
///

HWY_BEFORE_NAMESPACE();
namespace component_simd {
namespace HWY_NAMESPACE {

namespace hn = hwy::HWY_NAMESPACE;

// Hamming window coefficients (precomputed)
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

// DCT constants
static constexpr float CONST_SQRT2_2 = 0.707106781f;  // 1/sqrt(2)
static constexpr float CONST_C1 = 0.382683433f;       // cos(3*pi/8) - cos(pi/8)
static constexpr float CONST_C2 = 0.541196100f;       // cos(pi/8)
static constexpr float CONST_C3 = 1.306562965f;       // cos(pi/8) + cos(3*pi/8)

// 8x8 transpose using Highway - simple array-based approach
// This works correctly across all targets
template<class D, class V = hn::Vec<D>>
HWY_INLINE void Transpose8x8(D d, V &row0, V &row1, V &row2, V &row3,
                               V &row4, V &row5, V &row6, V &row7) {
  // Store to temporary array
  alignas(64) float temp[8][8];
  hn::Store(row0, d, temp[0]);
  hn::Store(row1, d, temp[1]);
  hn::Store(row2, d, temp[2]);
  hn::Store(row3, d, temp[3]);
  hn::Store(row4, d, temp[4]);
  hn::Store(row5, d, temp[5]);
  hn::Store(row6, d, temp[6]);
  hn::Store(row7, d, temp[7]);
  
  // Transpose in memory
  alignas(64) float transposed[8][8];
  for (int i = 0; i < 8; i++) {
    for (int j = 0; j < 8; j++) {
      transposed[i][j] = temp[j][i];
    }
  }
  
  // Load back
  row0 = hn::Load(d, transposed[0]);
  row1 = hn::Load(d, transposed[1]);
  row2 = hn::Load(d, transposed[2]);
  row3 = hn::Load(d, transposed[3]);
  row4 = hn::Load(d, transposed[4]);
  row5 = hn::Load(d, transposed[5]);
  row6 = hn::Load(d, transposed[6]);
  row7 = hn::Load(d, transposed[7]);
}

// Process one 8x8 DCT block using Highway SIMD
// Vectorizes across the 8 columns of the block
void ProcessBlock_HWY(class Line *inputLines[8], class Line *outputLines[8][8], int i) {
  // Use CappedTag to get at most 8 lanes
  const hn::CappedTag<float, 8> d;
  
  // Check if we have enough lanes for 8x8 vectorization
  if (hn::Lanes(d) < 8) {
    // Fallback to scalar for targets with < 8 lanes (e.g., SSE4 has 4)
    alignas(64) float block[8][8];
    for (int y = 0; y < 8; y++) {
      const float *src = inputLines[y]->Origin() + i;
      for (int x = 0; x < 8; x++) {
        block[y][x] = src[x];
      }
    }
    
    float sum = 0.0f;
    for (int y = 0; y < 8; y++) {
      for (int x = 0; x < 8; x++) {
        sum += block[y][x];
      }
    }
    float avg = sum / 64.0f;
    
    for (int y = 0; y < 8; y++) {
      for (int x = 0; x < 8; x++) {
        float win = window_table[y][x];
        block[y][x] = block[y][x] * win + avg * (1.0f - win);
      }
    }
    
    alignas(64) float temp[8][8];
    for (int y = 0; y < 8; y++) {
      float tmp0 = block[y][0] + block[y][7];
      float tmp7 = block[y][0] - block[y][7];
      float tmp1 = block[y][1] + block[y][6];
      float tmp6 = block[y][1] - block[y][6];
      float tmp2 = block[y][2] + block[y][5];
      float tmp5 = block[y][2] - block[y][5];
      float tmp3 = block[y][3] + block[y][4];
      float tmp4 = block[y][3] - block[y][4];
      
      float tmp10 = tmp0 + tmp3;
      float tmp13 = tmp0 - tmp3;
      float tmp11 = tmp1 + tmp2;
      float tmp12 = tmp1 - tmp2;
      
      temp[y][0] = tmp10 + tmp11;
      temp[y][4] = tmp10 - tmp11;
      
      float z1 = (tmp12 + tmp13) * CONST_SQRT2_2;
      temp[y][2] = tmp13 + z1;
      temp[y][6] = tmp13 - z1;
      
      tmp10 = tmp4 + tmp5;
      tmp11 = tmp5 + tmp6;
      tmp12 = tmp6 + tmp7;
      
      float z5 = (tmp10 - tmp12) * CONST_C1;
      float z2 = CONST_C2 * tmp10 + z5;
      float z4 = CONST_C3 * tmp12 + z5;
      float z3 = tmp11 * CONST_SQRT2_2;
      
      float z11 = tmp7 + z3;
      float z13 = tmp7 - z3;
      
      temp[y][5] = z13 + z2;
      temp[y][3] = z13 - z2;
      temp[y][1] = z11 + z4;
      temp[y][7] = z11 - z4;
    }
    
    alignas(64) float output[8][8];
    for (int x = 0; x < 8; x++) {
      float tmp0 = temp[0][x] + temp[7][x];
      float tmp7 = temp[0][x] - temp[7][x];
      float tmp1 = temp[1][x] + temp[6][x];
      float tmp6 = temp[1][x] - temp[6][x];
      float tmp2 = temp[2][x] + temp[5][x];
      float tmp5 = temp[2][x] - temp[5][x];
      float tmp3 = temp[3][x] + temp[4][x];
      float tmp4 = temp[3][x] - temp[4][x];
      
      float tmp10 = tmp0 + tmp3;
      float tmp13 = tmp0 - tmp3;
      float tmp11 = tmp1 + tmp2;
      float tmp12 = tmp1 - tmp2;
      
      output[0][x] = tmp10 + tmp11;
      output[4][x] = tmp10 - tmp11;
      
      float z1 = (tmp12 + tmp13) * CONST_SQRT2_2;
      output[2][x] = tmp13 + z1;
      output[6][x] = tmp13 - z1;
      
      tmp10 = tmp4 + tmp5;
      tmp11 = tmp5 + tmp6;
      tmp12 = tmp6 + tmp7;
      
      float z5 = (tmp10 - tmp12) * CONST_C1;
      float z2 = CONST_C2 * tmp10 + z5;
      float z4 = CONST_C3 * tmp12 + z5;
      float z3 = tmp11 * CONST_SQRT2_2;
      
      float z11 = tmp7 + z3;
      float z13 = tmp7 - z3;
      
      output[5][x] = z13 + z2;
      output[3][x] = z13 - z2;
      output[1][x] = z11 + z4;
      output[7][x] = z11 - z4;
    }
    
    for (int x = 0; x < 8; x++) {
      for (int y = 0; y < 8; y++) {
        *(outputLines[x][y]->Origin() + i) = output[x][y];
      }
    }
    return;
  }
  
  // SIMD path: Load one 8x8 block with each vector holding one row of 8 pixels
  using V = hn::Vec<decltype(d)>;
  V in[8];
  for (int y = 0; y < 8; y++) {
    in[y] = hn::LoadU(d, inputLines[y]->Origin() + i);
  }
  
  // Compute average of all 64 pixels
  V sum = hn::Add(hn::Add(hn::Add(in[0], in[1]), hn::Add(in[2], in[3])),
                  hn::Add(hn::Add(in[4], in[5]), hn::Add(in[6], in[7])));
  alignas(32) float sum_arr[8];
  hn::Store(sum, d, sum_arr);
  float avg = (sum_arr[0] + sum_arr[1] + sum_arr[2] + sum_arr[3] +
               sum_arr[4] + sum_arr[5] + sum_arr[6] + sum_arr[7]) / 64.0f;
  
  // Apply windowing using precomputed symmetric window vectors
  auto avg_vec = hn::Set(d, avg);
  auto one_vec = hn::Set(d, 1.0f);
  
  // Process symmetric rows: 0=7, 1=6, 2=5, 3=4
  for (int row_idx = 0; row_idx < 4; row_idx++) {
    auto win_vec = hn::LoadU(d, window_table[row_idx]);
    auto one_minus_win = hn::Sub(one_vec, win_vec);
    auto term = hn::Mul(avg_vec, one_minus_win);
    // in = in * win + avg * (1 - win)
    in[row_idx] = hn::MulAdd(in[row_idx], win_vec, term);
    in[7 - row_idx] = hn::MulAdd(in[7 - row_idx], win_vec, term);
  }
  
  // First DCT pass: process rows (vectorized across 8 columns)
  auto tmp0 = hn::Add(in[0], in[7]);
  auto tmp7 = hn::Sub(in[0], in[7]);
  auto tmp1 = hn::Add(in[1], in[6]);
  auto tmp6 = hn::Sub(in[1], in[6]);
  auto tmp2 = hn::Add(in[2], in[5]);
  auto tmp5 = hn::Sub(in[2], in[5]);
  auto tmp3 = hn::Add(in[3], in[4]);
  auto tmp4 = hn::Sub(in[3], in[4]);
  
  auto tmp10 = hn::Add(tmp0, tmp3);
  auto tmp13 = hn::Sub(tmp0, tmp3);
  auto tmp11 = hn::Add(tmp1, tmp2);
  auto tmp12 = hn::Sub(tmp1, tmp2);
  
  in[0] = hn::Add(tmp10, tmp11);
  in[4] = hn::Sub(tmp10, tmp11);
  
  auto z1 = hn::Mul(hn::Add(tmp12, tmp13), hn::Set(d, CONST_SQRT2_2));
  in[2] = hn::Add(tmp13, z1);
  in[6] = hn::Sub(tmp13, z1);
  
  tmp10 = hn::Add(tmp4, tmp5);
  tmp11 = hn::Add(tmp5, tmp6);
  tmp12 = hn::Add(tmp6, tmp7);
  
  auto z5 = hn::Mul(hn::Sub(tmp10, tmp12), hn::Set(d, CONST_C1));
  auto z2 = hn::MulAdd(hn::Set(d, CONST_C2), tmp10, z5);
  auto z4 = hn::MulAdd(hn::Set(d, CONST_C3), tmp12, z5);
  auto z3 = hn::Mul(tmp11, hn::Set(d, CONST_SQRT2_2));
  
  auto z11 = hn::Add(tmp7, z3);
  auto z13 = hn::Sub(tmp7, z3);
  
  in[5] = hn::Add(z13, z2);
  in[3] = hn::Sub(z13, z2);
  in[1] = hn::Add(z11, z4);
  in[7] = hn::Sub(z11, z4);
  
  // Transpose: convert from row-major to column-major for the second pass
  Transpose8x8(d, in[0], in[1], in[2], in[3], in[4], in[5], in[6], in[7]);
  
  // Second DCT pass (same operations)
  tmp0 = hn::Add(in[0], in[7]);
  tmp7 = hn::Sub(in[0], in[7]);
  tmp1 = hn::Add(in[1], in[6]);
  tmp6 = hn::Sub(in[1], in[6]);
  tmp2 = hn::Add(in[2], in[5]);
  tmp5 = hn::Sub(in[2], in[5]);
  tmp3 = hn::Add(in[3], in[4]);
  tmp4 = hn::Sub(in[3], in[4]);
  
  tmp10 = hn::Add(tmp0, tmp3);
  tmp13 = hn::Sub(tmp0, tmp3);
  tmp11 = hn::Add(tmp1, tmp2);
  tmp12 = hn::Sub(tmp1, tmp2);
  
  in[0] = hn::Add(tmp10, tmp11);
  in[4] = hn::Sub(tmp10, tmp11);
  
  z1 = hn::Mul(hn::Add(tmp12, tmp13), hn::Set(d, CONST_SQRT2_2));
  in[2] = hn::Add(tmp13, z1);
  in[6] = hn::Sub(tmp13, z1);
  
  tmp10 = hn::Add(tmp4, tmp5);
  tmp11 = hn::Add(tmp5, tmp6);
  tmp12 = hn::Add(tmp6, tmp7);
  
  z5 = hn::Mul(hn::Sub(tmp10, tmp12), hn::Set(d, CONST_C1));
  z2 = hn::MulAdd(hn::Set(d, CONST_C2), tmp10, z5);
  z4 = hn::MulAdd(hn::Set(d, CONST_C3), tmp12, z5);
  z3 = hn::Mul(tmp11, hn::Set(d, CONST_SQRT2_2));
  
  z11 = hn::Add(tmp7, z3);
  z13 = hn::Sub(tmp7, z3);
  
  in[5] = hn::Add(z13, z2);
  in[3] = hn::Sub(z13, z2);
  in[1] = hn::Add(z11, z4);
  in[7] = hn::Sub(z11, z4);
  
  // Store results to output lines
  for (int x = 0; x < 8; x++) {
    alignas(32) float temp[8];
    hn::Store(in[x], d, temp);
    
    // Unroll stores for better instruction scheduling
    *(outputLines[x][0]->Origin() + i) = temp[0];
    *(outputLines[x][1]->Origin() + i) = temp[1];
    *(outputLines[x][2]->Origin() + i) = temp[2];
    *(outputLines[x][3]->Origin() + i) = temp[3];
    *(outputLines[x][4]->Origin() + i) = temp[4];
    *(outputLines[x][5]->Origin() + i) = temp[5];
    *(outputLines[x][6]->Origin() + i) = temp[6];
    *(outputLines[x][7]->Origin() + i) = temp[7];
  }
}

} // namespace HWY_NAMESPACE
} // namespace component_simd
HWY_AFTER_NAMESPACE();

#if HWY_ONCE
namespace component_simd {

// Export function for dynamic dispatch
HWY_EXPORT(ProcessBlock_HWY);

// Wrapper function
void ProcessBlock(class Line *inputLines[8], class Line *outputLines[8][8], int i) {
  HWY_DYNAMIC_DISPATCH(ProcessBlock_HWY)(inputLines, outputLines, i);
}

} // namespace component_simd

/// Component::Component
Component::Component(void)
  : m_pOutputBuffer(nullptr), m_ulOutputStride(0)
  , m_ulY(0)
  , m_iCores(1)
  , m_pOriginal(nullptr)
{
  memset(m_pInputLines,0,sizeof(m_pInputLines));
  memset(m_pOutputLines,0,sizeof(m_pOutputLines));

  // Window function is precomputed in the Highway SIMD namespace (window_table)
  // No need to compute m_Window here - it was never used in the optimized code
}
///

/// Component::~Component
Component::~Component(void)
{
  int x,y;

  // Delete Line objects (they don't own their data)
  for(y = 0;y < 8;y++) {
    for(x = 0;x < 8;x++) {
      delete m_pOutputLines[x][y];
      m_pOutputLines[x][y] = nullptr;
    }
    delete m_pInputLines[y];
    m_pInputLines[y] = nullptr;
  }
  
  // Free the contiguous buffer
  if (m_pOutputBuffer) {
    free(m_pOutputBuffer);
    m_pOutputBuffer = nullptr;
  }
  
  delete m_pOriginal;
}
///

/// Component::AllocateBuffers
// Eagerly allocate all input and output buffers once width is known.
// This is called BEFORE any processing starts (before warmup loop).
// After this, all pointers are valid and never nullptr.
void Component::AllocateBuffers(uint32_t width)
{
  int x, y;
  uint32_t w = width - 8;  // Output DCT line width (8-pixel overlap removed)
  
  assert(width > 8);  // Sanity check
  assert(m_pInputLines[0] == nullptr);  // Should only be called once
  
  // Allocate all 8 input lines upfront
  for (y = 0; y < 8; y++) {
    m_pInputLines[y] = new Line(width);
  }
  
  // Allocate all 64 output lines in one contiguous buffer for cache efficiency
  // Calculate stride: line length rounded up to 64-byte (cache line) alignment
  size_t bytesPerLine = w * sizeof(float);
  m_ulOutputStride = ((bytesPerLine + 63) / 64) * 64 / sizeof(float);
  
  // Allocate single contiguous buffer for all 64 lines
  size_t totalBytes = m_ulOutputStride * sizeof(float) * 64;
  m_pOutputBuffer = static_cast<float*>(aligned_alloc(64, totalBytes));
  
  if (!m_pOutputBuffer) {
    // Fallback to regular new if aligned_alloc fails
    m_pOutputBuffer = new float[m_ulOutputStride * 64];
  }
  
  // Create Line wrapper objects that point into the contiguous buffer
  for (y = 0; y < 8; y++) {
    for (x = 0; x < 8; x++) {
      int index = x * 8 + y;
      float *lineStart = m_pOutputBuffer + (index * m_ulOutputStride);
      m_pOutputLines[x][y] = new Line(lineStart, w);
    }
  }
  
  // Allocate original line (always allocate for consistency)
  m_pOriginal = new Line(w);
}
///

/// Component::PushLine
// Push a new line into the DCT conversion. Note that this is an overcomplete
// representation, i.e. each line triggers a new DCT transformation.
// PRECONDITION: AllocateBuffers() must have been called first.
void Component::PushLine(Line &line)
{
  assert(line.LengthOf() > 8);
  assert(m_pInputLines[0] != nullptr);  // AllocateBuffers must be called first
  
  // Swap line into pre-allocated buffer at current position (no allocation needed)
  m_pInputLines[m_ulY]->Swap(line);
  
  // After 8 lines buffered (m_ulY == 7), run DCT and rotate
  if (m_ulY == 7) {
    class Line *top;
    int y;
    
    // Copy center line to m_pOriginal if needed for special weighting modes
#if defined(WEIGHT_MSE) || defined(WEIGHT_DELTA_E)
    uint32_t w = m_pOriginal->LengthOf();
    float *orig_ptr = m_pOriginal->Origin();
    const float *src_ptr = m_pInputLines[4]->Origin() + 4;
    memcpy(orig_ptr, src_ptr, w * sizeof(float));
#endif
    
    // Run DCT transformation (fills all 64 pre-allocated output bands)
    DCT();
    
    // Rotate input buffer (circular buffer, reuse oldest line slot)
    top = m_pInputLines[0];
    for (y = 0; y < 7; y++) {
      m_pInputLines[y] = m_pInputLines[y+1];
    }
    m_pInputLines[7] = top;
    // m_ulY stays at 7 (always triggers DCT from now on)
  } else {
    // Warmup phase: increment position until buffer is full
    m_ulY++;
  }
}
///

///Component::Run 
// Highway SIMD version - portable across x86 (SSE4/AVX2/AVX-512), ARM (NEON), etc.
void Component::Run(int offset,int mod)
{  
  int i,w; 
  w = m_pInputLines[0]->LengthOf();

  for(i = offset;i < w - 8;i+=mod) {
    component_simd::ProcessBlock(m_pInputLines, m_pOutputLines, i);
  }
}


void Component::DCT(void)
{
  // Run the "real" work now.
  SplitWork(m_iCores);
}
///

#endif // HWY_ONCE

