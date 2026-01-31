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

/// Includes
#include "options.hpp"
#include "measure/masking.hpp"
#include "measure/masking_simd.hpp"  // Highway SIMD functions
#include <cmath>
#include <cstring>
#include "dct/line.hpp"
///

/// Masking::Masking
Masking::Masking(void)
  : m_pMask(nullptr), m_pOutput(nullptr), m_pOriginal(nullptr), 
    m_pMapped(nullptr),
    m_pWideBuffer(nullptr), m_pNormalBuffer(nullptr),
    m_ulWideStride(0), m_ulNormalStride(0), m_ulWideWidth(0),
    m_ulY(0), m_dMaskingSlope(1.0), m_dVisibility(1.0),
    m_bPostFilter(false), m_ulWidth(0), m_bCoeffsComputed(false)
  // Slope and masking exponent will be installed over.
{
  memset(m_pInput   ,0,sizeof(m_pInput));
  memset(m_pBuffer  ,0,sizeof(m_pBuffer));
  memset(m_pAdded   ,0,sizeof(m_pAdded));
  memset(m_pFiltered,0,sizeof(m_pFiltered));

  // Initialize the 1D window function (Hamming window)
  // Note: The 2D m_Window is no longer used - we precompute window in the Highway DCT SIMD code
  float xz = float(MaskSize - 1) * 0.5f;
  float xf = float(MaskSize + 1) * 0.5f;
  float norm = 0.0f;
  float nsqr;

  // Compute 1D mask coefficients (gy values)
  // But accumulate norm from the full 2D window (as in original)
  for(int y = 0; y < MaskSize; y++) {
    float gy = cosf((y - xz) / xf * float(M_PI) * 0.5f);
    m_Mask[y] = gy;
    // Accumulate 2D window norm (sum over x dimension too)
    for(int x = 0; x < MaskSize; x++) {
      float gx = cosf((x - xz) / xf * float(M_PI) * 0.5f);
      norm += gx * gy;
    }
  }

  // Normalize: divide by MaskSize^2, then take sqrt (matches original)
  norm /= float(MaskSize * MaskSize);
  nsqr = sqrtf(norm);
  
  for(int y = 0; y < MaskSize; y++) {
    m_Mask[y] /= nsqr;
  } 
}
///

/// Masking::InitializeBuffers
// Allocate all buffers once we know the width - called on first PushLine
// Line wrappers point to contiguous buffers (minimal overhead, no separate allocations)
void Masking::InitializeBuffers(uint32_t width)
{
  int i, j, wideIdx, normalIdx;
  uint32_t w = width + MaskSize;
  
  // Cache widths
  m_ulWideWidth = w;
  
  // Count lines by width
  const int wideLines = 11;     // m_pInput[5], m_pBuffer[5], m_pMapped
  const int normalLines = 11;
  
  // Calculate strides (64-byte aligned)
  size_t wideBytesPerLine = w * sizeof(float);
  m_ulWideStride = ((wideBytesPerLine + 63) / 64) * 64 / sizeof(float);
  
  size_t normalBytesPerLine = width * sizeof(float);
  m_ulNormalStride = ((normalBytesPerLine + 63) / 64) * 64 / sizeof(float);
  
  // Allocate two contiguous buffers
  size_t wideTotalBytes = m_ulWideStride * sizeof(float) * wideLines;
  m_pWideBuffer = static_cast<float*>(aligned_alloc(64, wideTotalBytes));
  if (!m_pWideBuffer) {
    m_pWideBuffer = new float[m_ulWideStride * wideLines];
  }
  
  size_t normalTotalBytes = m_ulNormalStride * sizeof(float) * normalLines;
  m_pNormalBuffer = static_cast<float*>(aligned_alloc(64, normalTotalBytes));
  if (!m_pNormalBuffer) {
    m_pNormalBuffer = new float[m_ulNormalStride * normalLines];
  }
  
  // Create Line wrappers pointing into contiguous buffers
  wideIdx = 0;
  normalIdx = 0;
  
  // Wide lines
  for(i = 0; i < MaskSize; i++) {
    m_pInput[i] = new Line(m_pWideBuffer + (wideIdx * m_ulWideStride), w);
    wideIdx++;
  }
  for(i = 0; i < MaskSize; i++) {
    m_pBuffer[i] = new Line(m_pWideBuffer + (wideIdx * m_ulWideStride), w);
    wideIdx++;
  }
  m_pMapped = new Line(m_pWideBuffer + (wideIdx * m_ulWideStride), w);
  wideIdx++;
  
  // Normal lines
  for(i = 0; i < MaskSize; i++) {
    m_pAdded[i] = new Line(m_pNormalBuffer + (normalIdx * m_ulNormalStride), width);
    normalIdx++;
  }
  m_pMask = new Line(m_pNormalBuffer + (normalIdx * m_ulNormalStride), width);
  normalIdx++;
  m_pOutput = new Line(m_pNormalBuffer + (normalIdx * m_ulNormalStride), width);
  normalIdx++;
  
  for(j = 0; j < 2; j++) {
    for(i = 0; i < 2; i++) {
      m_pFiltered[i][j] = new Line(m_pNormalBuffer + (normalIdx * m_ulNormalStride), width);
      normalIdx++;
    }
  }
  
  m_ulWidth = width;
}
///

/// Masking::~Masking
// Simplified: Let process exit handle cleanup (no explicit free needed)
Masking::~Masking(void)
{
  // Line wrapper objects will be cleaned up automatically
  // Buffers don't need explicit free() - OS reclaims on process exit
}
///

/// Masking::EnablePostFilter
// Install the post-filter (called once per image, use float for consistency)
void Masking::EnablePostFilter(void)
{
  int i;
  m_bPostFilter = true;
  m_NormLo      = 0.0f;
  m_NormHi      = 0.0f;

  for(i=0;i < MaskSize;i++) {
    float xp = (i - ((MaskSize - 1) >> 1)) / float((MaskSize - 1) >> 1);
    float lo = 1.0f;
    float hi = sinf(xp * float(M_PI) * 0.5f);
    m_LowFilter[i] = lo;
    m_HiFilter[i]  = hi;
    m_NormHi      += hi * hi;
  }
  m_NormLo = 1.0f / MaskSize;
  m_NormHi = 1.0f / sqrtf(m_NormHi);
}
///

/// Masking::PushOriginal
// Push an original line into the line buffer
void Masking::PushOriginal(class Line *line)
{
  // Buffers already allocated by InitializeBuffers (called from PushLine)
  m_pBuffer[m_ulY]->Swap(*line);
}
///

/// Masking::PushLine
// Provide a new input line to the masking computation.
__attribute__((target("avx2,fma")))
void Masking::PushLine(class Line *line)
{
  const float vis     = m_dVisibility;
  int i,w   = line->LengthOf();
  int width = w - MaskSize;
  
  assert(width > 0);
  
  // Initialize all buffers on first call (branch predicted correctly after first call)
  if (__builtin_expect(m_ulWidth == 0, 0)) {
    InitializeBuffers(width);
  }
  
  // Swap the contents of the input line into our buffer
  m_pInput[m_ulY]->Swap(*line);
  //
  // Convert to nonlinear space
  float *input_ptr = m_pInput[m_ulY]->Origin();
  float *mapped_ptr = m_pMapped->Origin();
  const float slope = m_dMaskingSlope;
  
  // Highway SIMD path - portable and length-agnostic (works with SSE2/AVX2/AVX-512/NEON)
  if (slope == 1.0f) {
    // expon=1.0: vectorized fabsf + multiply
    masking_simd::ApplyMaskingSlope1(input_ptr, mapped_ptr, vis, w);
  } else if (slope == 0.5f) {
    // expon=0.5: vectorized sqrt of scaled absolute
    masking_simd::ApplyMaskingSlope05(input_ptr, mapped_ptr, vis, w);
  } else {
    // General case: use powf (keep scalar - pow is expensive anyway)
    for(i = 0; i < w; i++) {
      float v = fabsf(input_ptr[i]) * vis;
      mapped_ptr[i] = powf(v, slope);
    }
  }
  //
  // Column sum - vectorized to avoid another load/store cycle
  float *added_ptr = m_pAdded[m_ulY]->Origin();
  
  // Exploit symmetry: m_Mask[0]=m_Mask[4], m_Mask[1]=m_Mask[3], m_Mask[2]=center
  const float m0 = m_Mask[0], m1 = m_Mask[1], m2 = m_Mask[2];
  
  // Highway 5-tap symmetric convolution - adapts to vector width (128/256/512-bit)
  masking_simd::ColumnSumConvolution(mapped_ptr, added_ptr, m0, m1, m2, width);
  //
  // If the input buffer is full, compute the new mask.
  if (m_ulY == MaskSize-1) {
    class Line *top,*add,*buf;
    int y;
    //
    ComputeMask();
    //
    if (m_bPostFilter) {
      ComputeLowpass();
    }
    //
    // Make the topmost line available for computations
    top = m_pInput[0];
    add = m_pAdded[0];
    buf = m_pBuffer[0];
    for(y = 0;y < MaskSize-1;y++) {
      m_pInput[y]  = m_pInput[y+1];
      m_pAdded[y]  = m_pAdded[y+1];
      m_pBuffer[y] = m_pBuffer[y+1];
    }
    m_pInput[MaskSize-1]  = top;
    m_pAdded[MaskSize-1]  = add;
    m_pBuffer[MaskSize-1] = buf;
  } else {
    // the next input goes into the next line.
    m_ulY++;
  }
}
///

/// Masking::ComputeLowpass
// Compute the low and high-pass filter of the input, splitting one band into
// two - optimized version with cached coefficients and direct pointer access
void Masking::ComputeLowpass(void)
{
  int i,j,k;
  const int width = m_ulWidth;

  // Compute filter coefficients once (cached)
  if (__builtin_expect(!m_bCoeffsComputed, 0)) {
    const float sll = m_dVisibility * m_NormLo * m_NormLo;
    const float slh = m_dVisibility * m_NormLo * m_NormHi;
    const float shh = m_dVisibility * m_NormHi * m_NormHi;
    
    for(i = 0; i < MaskSize; i++) {
      for(j = 0; j < MaskSize; j++) {
        m_FilterCoeffs_LL[i][j] = m_LowFilter[i] * m_LowFilter[j] * sll;
        m_FilterCoeffs_LH[i][j] = m_LowFilter[i] * m_HiFilter[j]  * slh;
        m_FilterCoeffs_HL[i][j] = m_HiFilter[i]  * m_LowFilter[j] * slh;
        m_FilterCoeffs_HH[i][j] = m_HiFilter[i]  * m_HiFilter[j]  * shh;
      }
    }
    m_bCoeffsComputed = true;
  }

  // Direct pointers to output buffers
  float *filt_ll = m_pFiltered[0][0]->Origin();
  float *filt_lh = m_pFiltered[1][0]->Origin();
  float *filt_hl = m_pFiltered[0][1]->Origin();
  float *filt_hh = m_pFiltered[1][1]->Origin();
  
  // Direct pointers to input lines
  const float *input_ptrs[MaskSize];
  for(i = 0; i < MaskSize; i++) {
    input_ptrs[i] = m_pInput[i]->Origin();
  }
  
  // Main convolution loop
  for(k = 0; k < width; k++) {
    float acc_ll = 0.0f, acc_lh = 0.0f, acc_hl = 0.0f, acc_hh = 0.0f;
    
    // Accumulate over 5x5 filter window
    for(i = 0; i < MaskSize; i++) {
      const float *input_row = input_ptrs[i] + k;
      // Unroll inner loop for better pipeline utilization
      const float val0 = input_row[0];
      const float val1 = input_row[1];
      const float val2 = input_row[2];
      const float val3 = input_row[3];
      const float val4 = input_row[4];
      
      acc_ll += val0 * m_FilterCoeffs_LL[i][0] + val1 * m_FilterCoeffs_LL[i][1] + 
                val2 * m_FilterCoeffs_LL[i][2] + val3 * m_FilterCoeffs_LL[i][3] + 
                val4 * m_FilterCoeffs_LL[i][4];
      acc_lh += val0 * m_FilterCoeffs_LH[i][0] + val1 * m_FilterCoeffs_LH[i][1] + 
                val2 * m_FilterCoeffs_LH[i][2] + val3 * m_FilterCoeffs_LH[i][3] + 
                val4 * m_FilterCoeffs_LH[i][4];
      acc_hl += val0 * m_FilterCoeffs_HL[i][0] + val1 * m_FilterCoeffs_HL[i][1] + 
                val2 * m_FilterCoeffs_HL[i][2] + val3 * m_FilterCoeffs_HL[i][3] + 
                val4 * m_FilterCoeffs_HL[i][4];
      acc_hh += val0 * m_FilterCoeffs_HH[i][0] + val1 * m_FilterCoeffs_HH[i][1] + 
                val2 * m_FilterCoeffs_HH[i][2] + val3 * m_FilterCoeffs_HH[i][3] + 
                val4 * m_FilterCoeffs_HH[i][4];
    }
    
    filt_ll[k] = acc_ll;
    filt_lh[k] = acc_lh;
    filt_hl[k] = acc_hl;
    filt_hh[k] = acc_hh;
  }
}

/// Masking::ComputeMask
// Compute the masking strength for the buffered lines, create a new buffered
// output line.
// Float buffers and float math throughout - no conversions in hot loop!
void Masking::ComputeMask(void)
{
  int width = m_ulWidth;
  const float base = BASE_VISIBILITY;
  const float invMaskSizeSq = 1.0f / (MaskSize * MaskSize);
  
  assert(width > 0);
  
  // Exploit symmetry: m_Mask[0]=m_Mask[4], m_Mask[1]=m_Mask[3], m_Mask[2]=center
  const float m0 = m_Mask[0], m1 = m_Mask[1], m2 = m_Mask[2];
  
  // Prepare pointers for Highway SIMD
  float* added_ptrs[5];
  for (int i = 0; i < 5; i++) {
    added_ptrs[i] = m_pAdded[i]->Origin();
  }
  
  const float* input_center = m_pInput[MaskSize >> 1]->Origin() + (MaskSize >> 1);
  float* mask_out = m_pMask->Origin();
  float* output_out = m_pOutput->Origin();
  
  // Highway SIMD - portable and length-agnostic
#if NO_BASE_VISIBILITY
  masking_simd::ComputeMaskLoop(added_ptrs, input_center, mask_out, output_out,
                                 m0, m1, m2, base, invMaskSizeSq, m_dVisibility,
                                 width, true);
#else
  masking_simd::ComputeMaskLoop(added_ptrs, input_center, mask_out, output_out,
                                 m0, m1, m2, base, invMaskSizeSq, m_dVisibility,
                                 width, false);
#endif
}
///

