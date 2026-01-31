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
#include "measure/pooling.hpp"
#include "img/image.hpp"
#include "img/imgwriter.hpp"
#include "measure/masking.hpp"
#include "dct/line.hpp"
#include <cstdio>
#include <cmath>
#include "measure/ediff.hpp"

// Highway for portable SIMD
#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "measure/pooling.cpp"
#include <hwy/foreach_target.h>
#include <hwy/highway.h>

HWY_BEFORE_NAMESPACE();
namespace pooling_simd {
namespace HWY_NAMESPACE {

namespace hn = hwy::HWY_NAMESPACE;

// Highway SIMD version of MeasureInBand
void MeasureInBandSIMD_HWY(float *ref_ptr, float *dst_ptr,
                            const float *refmask_ptr, const float *dstmask_ptr,
                            float *error_ptr, int simd_limit, float visbase) {
  const hn::ScalableTag<float> d;
  const size_t N = hn::Lanes(d);
  
  auto vvisbase = hn::Set(d, visbase);
  auto vzero = hn::Zero(d);
  
  // Process SIMD chunks
  for (int i = 0; i < simd_limit; i += N) {
    // Load masks
    auto vrefmask = hn::LoadU(d, refmask_ptr + i);
    auto vdstmask = hn::LoadU(d, dstmask_ptr + i);
    
    // vis = max(refmask, dstmask) * visbase
    auto vvis = hn::Max(vrefmask, vdstmask);
    vvis = hn::Mul(vvis, vvisbase);
    
    // err = refline - dstline
    auto vref = hn::LoadU(d, ref_ptr + i);
    auto vdst = hn::LoadU(d, dst_ptr + i);
    auto verr = hn::Sub(vref, vdst);
    
    // err = fabs(err) * vis
    verr = hn::Abs(verr);
    verr = hn::Mul(verr, vvis);
    
    // Compute err = exp(-pow(err, 3.5))
    // pow(err, 3.5) = err^2 * err * sqrt(err)
    auto verr2 = hn::Mul(verr, verr);      // err^2
    auto verr3 = hn::Mul(verr2, verr);     // err^3
    auto vsqrt_err = hn::Sqrt(verr);       // sqrt(err)
    auto verr_pow = hn::Mul(verr3, vsqrt_err); // err^3.5
    
    // Negate for exp(-x)
    verr_pow = hn::Neg(verr_pow);
    
    // Create mask for err > 0
    auto vmask = hn::Gt(verr, vzero);
    
    // Compute exp() - must do scalar for now (no SIMD exp in Highway)
    alignas(64) float err_pow_array[hn::MaxLanes(d)];
    alignas(64) float exp_result[hn::MaxLanes(d)];
    hn::StoreU(verr_pow, d, err_pow_array);
    
    // Compute exp for each element
    for (size_t j = 0; j < N; j++) {
      exp_result[j] = expf(err_pow_array[j]);
    }
    
    auto vexp_result = hn::LoadU(d, exp_result);
    
    // errorline[i] *= exp_result (where err > 0, else keep unchanged)
    auto verror = hn::LoadU(d, error_ptr + i);
    auto vresult = hn::Mul(verror, vexp_result);
    
    // Blend: use vresult where mask is true, verror where false
    vresult = hn::IfThenElse(vmask, vresult, verror);
    
    hn::StoreU(vresult, d, error_ptr + i);
  }
}

} // namespace HWY_NAMESPACE
} // namespace pooling_simd
HWY_AFTER_NAMESPACE();

#if HWY_ONCE
namespace pooling_simd {

// Export function for dynamic dispatch
HWY_EXPORT(MeasureInBandSIMD_HWY);

// Wrapper function
void MeasureInBandSIMD(float *ref_ptr, float *dst_ptr,
                       const float *refmask_ptr, const float *dstmask_ptr,
                       float *error_ptr, int simd_limit, float visbase) {
  HWY_DYNAMIC_DISPATCH(MeasureInBandSIMD_HWY)(ref_ptr, dst_ptr, refmask_ptr, dstmask_ptr,
                                               error_ptr, simd_limit, visbase);
}

} // namespace pooling_simd

///

/// Quantization tables
#ifdef AHUMADA
/*
 Perceptual quantization matrices from:

A VISUAL DETECTION MODEL FOR DCT QUANTIZATION

Albert J. Ahumada, Jr. and Andrew B. Watson
NASA Ames Research Center
Moffett Field, California
and
Heidi A. Peterson
IBM T. J. Watson Research Center
Yorktown Heights, New York 
*/
static const unsigned int luma_tbl[] = {
  25,  11,  11,  12,  15,  19,  25,  32,
  11,  13,  10,  10,  12,  15,  19,  24,
  11,  10,  14,  14,  16,  18,  22,  27,
  12,  10,  14,  18,  21,  24,  28,  33,
  15,  12,  16,  21,  26,  31,  36,  42,
  19,  15,  18,  24,  31,  38,  45,  53,
  25,  19,  22,  28,  36,  45,  55,  65,
  32,  24,  27,  33,  42,  53,  65,  77
};

static const unsigned int cr_tbl[] = {
  21,  21,  41,  45,  55,  71,  92, 120,
  21,  37,  39,  38,  44,  55,  70,  89,
  41,  39,  51,  54,  59,  69,  83, 103,
  45,  38,  54,  69,  80,  91, 106, 126,
  55,  44,  59,  80, 100, 117, 136, 158,
  71,  55,  69,  91, 117, 144, 170, 198,
  92,  70,  83, 106, 136, 170, 206, 243,
 120,  89, 103, 126, 158, 198, 243, 290
};

static const unsigned int cb_tbl[] = {
  45,  43, 103, 114, 141, 181, 236, 306,
  43,  78,  99,  97, 113, 140, 178, 228,
 103,  99, 130, 138, 150, 175, 212, 262,
 114,  97, 138, 176, 203, 232, 270, 321,
 141, 113, 150, 203, 254, 299, 347, 403,
 181, 140, 175, 232, 299, 367, 434, 505,
 236, 178, 212, 270, 347, 434, 525, 619,
 306, 228, 262, 321, 403, 505, 619, 739
};
#else
static const unsigned int luma_tbl[] = {
  16,  11,  10,  16,  24,  40,  51,  61,
  12,  12,  14,  19,  26,  58,  60,  55,
  14,  13,  16,  24,  40,  57,  69,  56,
  14,  17,  22,  29,  51,  87,  80,  62,
  18,  22,  37,  56,  68, 109, 103,  77,
  24,  35,  55,  64,  81, 104, 113,  92,
  49,  64,  78,  87, 103, 121, 120, 101,
  72,  92,  95,  98, 112, 100, 103,  99
};

/*
** Note that the following table is designed
** for color subsampling, i.e. frequencies
** are only half of the non-subsampled
** image
*/
static const unsigned int sub_chroma_tbl[] = {
  17,  18,  24,  47,  99,  99,  99,  99,
  18,  21,  26,  66,  99,  99,  99,  99,
  24,  26,  56,  99,  99,  99,  99,  99,
  47,  66,  99,  99,  99,  99,  99,  99,
  99,  99,  99,  99,  99,  99,  99,  99,
  99,  99,  99,  99,  99,  99,  99,  99,
  99,  99,  99,  99,  99,  99,  99,  99,
  99,  99,  99,  99,  99,  99,  99,  99
};

static const unsigned int chroma_tbl[] = {
  17,  24,  47,  99,  99,  99,  99,  99,
  24,  26,  66,  99,  99,  99,  99,  99,
  47,  66,  99,  99,  99,  99,  99,  99,
  99,  99,  99,  99,  99,  99,  99,  99,
  99,  99,  99,  99,  99,  99,  99,  99,
  99,  99,  99,  99,  99,  99,  99,  99,
  99,  99,  99,  99,  99,  99,  99,  99,
  99,  99,  99,  99,  99,  99,  99,  99
}; 
#endif

//
// Output of the DCT must be further divided by the
// numbers below to normalize it correctly.
static const float dct_scale[] = {
  1.0f, 1.387039845f, 1.306562965f, 1.175875602f,
  1.0f, 0.785694958f, 0.541196100f, 0.275899379f
};

//
// Additional correction factors that make the transformation
// equal-frequency dependent. It isn't currently due to the window
// function. What we see here is the l^2 norm squared under the
// transformation on an input that is given by a constant-frequency
// DCT pattern.
static const float norms[8][8] = {
  {1.0f,0.742298f,1.22273f,1.24434f,1.2528f,1.26126f,1.28286f,1.76329f},
  {0.742298f,0.439822f,0.724487f,0.737285f,0.742298f,0.747311f,0.760109f,1.04477f},
  {1.22273f,0.724487f,1.19339f,1.21448f,1.22273f,1.23099f,1.25207f,1.72098f},
  {1.24434f,0.737285f,1.21448f,1.23593f,1.24434f,1.25274f,1.27419f,1.75139f},
  {1.2528f,0.742298f,1.22273f,1.24434f,1.2528f,1.26126f,1.28286f,1.76329f},
  {1.26126f,0.747311f,1.23099f,1.25274f,1.26126f,1.26977f,1.29152f,1.7752f},
  {1.28286f,0.760109f,1.25207f,1.27419f,1.28286f,1.29152f,1.31364f,1.8056f},
  {1.76329f,1.04477f,1.72098f,1.75139f,1.76329f,1.7752f,1.8056f,2.48181f},
};
///

/// Pooling::Pooling
Pooling::Pooling(void)
{
}
///

/// Pooling::Pooling
Pooling::~Pooling(void)
{
}
///

/// Pooling::Run
void Pooling::Run(int offset,int mod)
{   
  int cnt   = 0;
  int comps = m_iComponents;
  int logme = m_iLogging;
  int x,y,c;
  
  for(c = 0;c < comps;c++) {
    for(y = 0;y < 8;y++) {
      for(x = 0;x < 8;x++) {
	{
	  if (cnt % mod == offset) {
	    // After warmup, all DCT bands are guaranteed valid (never nullptr)
	    class Line *refline = m_pReference->GetDCTBand(x,y,c);
	    class Line *dstline = m_pDistorted->GetDCTBand(x,y,c);
	    
	    if (logme > 1) {
	      m_DCTImgs1[x][y][c].WriteLine(refline);
	      m_DCTImgs2[x][y][c].WriteLine(dstline);
	    }
	    m_Weights1[x][y][c].PushLine(refline);
	    m_Weights2[x][y][c].PushLine(dstline);
	  }
	  cnt++;
	}
      }
    }
  }
}
///

// Fast approximation of exp(x) for x in [-10, 0] (our range after -pow())
// Uses polynomial approximation - much faster than libm exp()
static inline float fast_exp_negative(float x) {
  // Clamp to reasonable range
  if (x < -10.0f) return 0.0f;
  if (x > 0.0f) return 1.0f;
  
  // Use exp(x) ≈ 2^(x/ln(2)) and compute 2^k quickly
  // For negative x, we can use: exp(x) ≈ 1 / (1 - x + x^2/2 - x^3/6 + x^4/24)
  // Simplified 4th order polynomial approximation
  const float c1 = 1.0f;
  const float c2 = 1.0f;
  const float c3 = 0.5f;
  const float c4 = 0.16666667f;
  const float c5 = 0.041666667f;
  
  float x2 = x * x;
  float x3 = x2 * x;
  float x4 = x2 * x2;
  
  return c1 + x * c2 + x2 * c3 + x3 * c4 + x4 * c5;
}

/// Pooling::MeasureInBand
// Measure the visibility in a band
// and adjust the error map.
// Highway SIMD-optimized version - portable across x86 (SSE4/AVX2/AVX-512), ARM (NEON), etc.
void Pooling::MeasureInBand(class Line *refline,class Line *dstline,
			    const class Line *refmask,const class Line *dstmask,
			    int w,int logme,float visbase,
			    class Line *errorline)
{ 
  int i;
  
  // Get direct pointers for faster access
  float *ref_ptr = refline->Origin();
  float *dst_ptr = dstline->Origin();
  const float *refmask_ptr = refmask->Origin();
  const float *dstmask_ptr = dstmask->Origin();
  float *error_ptr = errorline->Origin();
  
#ifdef AHUMADA
  // AHUMADA_FACTOR=1.0, DETECT_THREAS=1.0, AHUMADA_EXPONENT=3.5
  // So: err = exp(-pow(err, 3.5)) = exp(-err^3.5)
  // pow(x, 3.5) = x^3 * sqrt(x) = x^2 * x * sqrt(x)
#endif
  
  // SIMD loop: process vectors at a time
  // If logme > 2 (verbose diagnostic output), fall back to scalar path to avoid branching
  // Since buffers are 64-byte aligned with padding, the last chunk can safely extend
  // into the padding zone (padding values are never read back)
  int simd_limit = (logme > 2) ? 0 : w;
  if (simd_limit > 0) {
    pooling_simd::MeasureInBandSIMD(ref_ptr, dst_ptr, refmask_ptr, dstmask_ptr,
                                     error_ptr, simd_limit, visbase);
  }
  
  /* Scalar tail loop NO LONGER NEEDED due to padding (when logme <= 2)!
   * 
   * When logme <= 2, all Line buffers have 64-byte aligned padding, allowing
   * us to safely process full SIMD chunks beyond 'w'. The padding values are
   * garbage but never read back since all later code respects 'w' as the true width.
   * 
   * When logme > 2, simd_limit is 0 and we use the full scalar path below for
   * diagnostic output (which involves branches that would hurt SIMD performance).
   * 
   * Original scalar code (kept for logme > 2 case):
   */
  for(i = (logme > 2) ? 0 : simd_limit; i < w; i++) {
    float vis    = 1.0;
    float visref = 1.0;
    float visdst = 1.0;
    float err;
    //
    // Get the relative visibility of the coefficients due to masking.
    // The DC part is not masked ("threshold elevation")
    visref = refmask->Get(i);
    visdst = dstmask->Get(i);
    // Artefacts in the compressed image should not impact the visibility, and
    // reduced visibility due to artefacts in the original not available in the
    // compressed image are also visible. For that reason, set the visibility to
    // the max of the two, or the elevation to the min of the two.
    vis  = (visref > visdst)?(visref):(visdst);
    vis *= visbase;
    //
    // Get the total square error and compute the detectiong probability.
    // Note that the masking algorithm already multiplied refline and dstline
    // with the DCT output factor and the CSF factor.
    err = refline->Get(i) - dstline->Get(i);
    err = fabsf(err); // DCTPSNR uses here err*err, but the difference gets visible, not its square
    err = err * vis;
    assert(err >= 0.0f);
    if (err > 0.0f) {
      // The visibility threshold is normalized such that
      // 1.0 is approximately 50% detection probabililty.
      // However, 1.0 is the "noticability" threshold, thus
      // upscale.
#ifdef AHUMADA
      err = expf(AHUMADA_FACTOR * -powf(DETECT_THREAS * err,AHUMADA_EXPONENT)); // 1.0 - detection probability.
#else
      err = expf(-powf(2.0f * err,3.5f)); // 1.0 - detection probability.
#endif
      assert(err >= 0.0 && err <= 1.0);
      //
      // Include in the total error detection probability:
      // Probability of not detecting an error requires that no band
      // detects an error.
      errorline->At(i) *= err;
      if (logme > 2) {
	refline->At(i) -= dstline->Get(i);
      }
    }
  }
}
///

/// Pooling::Measure
// Measure the difference between the two images, return the result.
// Optionally create output log images.
float Pooling::Measure(class Image *ref ,class Image *dist,
			class Image *sref,class Image *sdist,
			int cores,int logme)
{
  int x,y,c,i;
  const int offset = 8 + Masking::MaskSize;
  int w = ref->WidthOf()  - offset;
  int h = ref->HeightOf() - offset;
  int l = 0;
  int comps = ref->ComponentCountOf();
  class Line errorline(w);
  class Line backref(w),distref(w);
  float error = 0.0;
  //
  //
  m_iComponents = comps;
  m_iLogging    = logme;
  m_pReference  = ref;
  m_pDistorted  = dist;
  //
  // If saliency maps are defined, read the first lines to correct for the offsets.
  if (sref && sdist) {
    for(int si = 0; si < (offset >> 1); si++) {
      sref->ReadNextPFMLine();
      sdist->ReadNextPFMLine();
    }
  }
  //
  // Create the output images for logging if this is enabled.
  if (logme) {
#ifdef WEIGHT_MSE
    m_ErrorMap.OpenPGM(w,h,1023.0,0.0,"errors.pgm");
#else
    m_ErrorMap.OpenPGM(w,h,255.0,0.0,"errors.pgm");
#endif
    for(c = 0;c < comps;c++) {
      for(y = 0;y < 8;y++) {
	for(x = 0;x < 8;x++) {
	  float img_offset, img_scale;
	  if (x == 0 && y == 0) {
	    img_offset = 0.0f; 
	    img_scale  = 255.0f/8.0f; // show the average luminance, not the nominal low-pass
	  } else {
	    img_offset = 128.0f;
	    img_scale  = 128.0f;
	  }
	  // Correct for the wrong output scaling of the exponent.
	  img_scale /= (8.0f * dct_scale[x] * dct_scale[y] * sqrtf(norms[x][y]));
	  m_DCTImgs1[x][y][c].OpenPGM(ref->WidthOf() - 8,ref->HeightOf() - 8,img_scale,img_offset,"dct1_%d_%d_%d.pgm",x,y,c);
	  m_DCTImgs2[x][y][c].OpenPGM(dist->WidthOf() - 8,dist->HeightOf() - 8,img_scale,img_offset,"dct2_%d_%d_%d.pgm",x,y,c);
	  m_DCTDiff[x][y][c].OpenPGM(w,h,128.0f,128.0f,"dctdiff_%d_%d_%d.pgm",x,y,c);
	  m_MaskImgs1[x][y][c].OpenPGM(w,h,255.0f,0.0f,"mask1_%d_%d_%d.pgm",x,y,c);
	  m_MaskImgs2[x][y][c].OpenPGM(w,h,255.0f,0.0f,"mask2_%d_%d_%d.pgm",x,y,c);
#ifdef EXTENDED_FILTER
	  if (x == 0 && y == 0 && c == 0) {
	    // DCT scale correction is already done in the visibility computation
	    // of the mask, do not repeat.
	    float scale = luma_tbl[0] / 255.0f;
	    m_LowPasses1[0][0].OpenPGM(w,h,scale * 255.0f / 8.0f,0.0f   ,"lowpass1_%d_%d.pgm",0,0);
	    m_LowPasses2[0][0].OpenPGM(w,h,scale * 255.0f / 8.0f,0.0f   ,"lowpass2_%d_%d.pgm",0,0);
	    m_LowPasses1[0][1].OpenPGM(w,h,scale * 128.0f      ,128.0f ,"lowpass1_%d_%d.pgm",0,1);
	    m_LowPasses2[0][1].OpenPGM(w,h,scale * 128.0f      ,128.0f ,"lowpass2_%d_%d.pgm",0,1);
	    m_LowPasses1[1][0].OpenPGM(w,h,scale * 128.0f      ,128.0f ,"lowpass1_%d_%d.pgm",1,0);
	    m_LowPasses2[1][0].OpenPGM(w,h,scale * 128.0f      ,128.0f ,"lowpass2_%d_%d.pgm",1,0);
	    m_LowPasses1[1][1].OpenPGM(w,h,scale * 128.0f      ,128.0f ,"lowpass1_%d_%d.pgm",1,1);
	    m_LowPasses2[1][1].OpenPGM(w,h,scale * 128.0f      ,128.0f ,"lowpass2_%d_%d.pgm",1,1);
	  }
#endif
	}
      }
      {
	uint32_t wo = ref->WidthOf()  - 8;
	uint32_t ho = ref->HeightOf() - 8;
	float img_offset = (c == 0) ? 0.0f : 128.0f;
	m_SrcImgs1[c].OpenPGM(wo,ho,255.0f,img_offset,"refimg_%d.pgm",c);
	m_SrcImgs2[c].OpenPGM(wo,ho,255.0f,img_offset,"dstimg_%d.pgm",c);
      }
    }
  }

  //
  // Setup the masking and visibility weights for the bands.
  for(c = 0;c < comps;c++) {
    for(y = 0;y < 8;y++) {
      for(x = 0;x < 8;x++) {
	float vis,expon;
	int xs = x; // << 1;
	int ys = y; // << 1;
	if (xs > 7) xs = 7;
	if (ys > 7) ys = 7;
	//
	// First integrate the scaling of the DCT. Input must be scaled by 8, then by the
	// norm correction factor.
	vis = 1.0f / (8.0f * dct_scale[x] * dct_scale[y] * sqrtf(norms[x][y]));
	//
	// The idea of the next scaling is to bring the visibility just to
	// the "just noticable difference", which is by an 8bpp input given
	// by the quantization coefficients of the JPEG compression. All other
	// coefficients are then scaled to the visibility threshold of the DC
	// luminance.
	vis *= 255.0f / 16.0f;
	//
	// Further, include the base visibility of the frequency band.
	switch(c) {
	case 0: // Y
	  vis *= 1.0f   * 16.0f / luma_tbl  [x + (y << 3)];
	  break;
	case 1: // Cb
#ifdef AHUMADA
	  vis *= 0.088f * 45.0f / cb_tbl    [xs + (ys << 3)];
#else
	  vis *= 0.088f * 16.0f / chroma_tbl[x + (y << 3)];
#endif
	  break;
	case 2: // Cr
#ifdef AHUMADA
	  vis *= 0.278f * 21.0f / cr_tbl    [xs + (ys << 3)];
#else
	  vis *= 0.278f * 16.0f / chroma_tbl[x + (y << 3)];
#endif
	  break;
	}
	//
	// Get the masking exponent. This is also frequency dependent.
	switch(x+y) {
	case 0:
	  expon = 0.7f;
	  break;
	case 1:
	case 2:
	  expon = 0.8f;
	  break;
	case 3:
	  expon = 0.9f;
	  break;
	default:
	  expon = 1.0;
	  break;
	}
	//
	// And install the settings.
	m_Weights1[x][y][c].SetVisibility(vis);
	m_Weights2[x][y][c].SetVisibility(vis);
	m_Weights1[x][y][c].SetMaskingExponent(expon);
	m_Weights2[x][y][c].SetMaskingExponent(expon);
      }
    }
  }
  //
#ifdef EXTENDED_FILTER
  // The luminance low-pass resolution is insufficient, enlarge it by a post-filter.
  m_Weights1[0][0][0].EnablePostFilter();
  m_Weights2[0][0][0].EnablePostFilter();
#endif
  //
  // EAGER ALLOCATION: Allocate all Component buffers upfront (before any computation).
  // This ensures all pointers are valid (never nullptr) after warmup completes.
  // Width is known after OpenPNM, so we can allocate everything now.
  ref->AllocateAllBuffers();
  dist->AllocateAllBuffers();
  //
  // Explicit warmup phase: Fill DCT and masking buffers before main processing.
  // DCT needs 8 lines, masking needs 5 more (total 13 lines until all outputs valid).
  // After warmup, all GetDCTBand() and GetMask() calls return valid pointers (never nullptr).
  const int DCT_WARMUP = 8;
  const int MASK_WARMUP = 5;
  const int TOTAL_WARMUP = DCT_WARMUP + MASK_WARMUP - 1;  // 12 lines (they overlap by 1)
  
  int warmup_count = 0;
  while (warmup_count < TOTAL_WARMUP && !ref->ImageDone()) {
    ref->ReadNextLine();
    dist->ReadNextLine();
    
    // After DCT warmup (8 lines), start running masking (which needs 5 lines to warm up)
    if (warmup_count >= DCT_WARMUP - 1) {
      SplitWork(cores);  // Pushes lines into masking
    }
    
    warmup_count++;
  }
  //
  // Now for the main loop. All buffers are warmed up, no nullptr checks needed.
  while(!ref->ImageDone()) {
    bool collected = false;
    ref->ReadNextLine();
    dist->ReadNextLine(); 
    //
    //
    // Reset the detection probability for this line.
    for(i = 0;i < w;i++) {
      errorline.At(i) = 1.0f;
    }
    //
#if defined(WEIGHT_MSE) || defined(WEIGHT_DELTA_E)
    for(c = 0;c < comps;c++) {
      class Line *line;
      if ((line = m_pReference->GetOriginal(c))) {
	m_Weights1[0][0][c].PushOriginal(line);
      }
      if ((line = m_pDistorted->GetOriginal(c))) {
	m_Weights2[0][0][c].PushOriginal(line);
      }
    }
#endif
    if (logme > 1) {
      for(c = 0;c < comps;c++) {
	class Line *line;
	if ((line = m_pReference->GetOriginal(c))) {
	  m_SrcImgs1[c].WriteLine(line);
	}
	if ((line = m_pDistorted->GetOriginal(c))) {
	  m_SrcImgs2[c].WriteLine(line);
	}
      }
    }
    //
    // First run the masking algorithm.
    SplitWork(cores);
    //
    // Write the DCT components out for logging purposes.
    for(c = 0;c < comps;c++) {
      // Re-generate the mask for the DC-band. The model is here that the neighbour-band(s) generate the
      // visual mask for this band as otherwise artifacts around boundaries become visible.
      {
	class Line *mskdest  = m_Weights1[0][0][0].GetMask();
	class Line *msk1line = m_Weights1[0][1][0].GetMask();
	class Line *msk2line = m_Weights1[1][0][0].GetMask();
	class Line *msk3line = m_Weights1[1][1][0].GetMask();
	//
	if (mskdest) {
	  for(i = 0;i < w;i++) {
	    mskdest->At(i) = msk1line->Get(i) * 0.5f + msk2line->Get(i) * 0.5f + msk3line->Get(i) * 0.25f;
	  }
	}
      }
      //
      {
	// After warmup, all masks guaranteed valid (never nullptr)
	class Line *mskdest  = m_Weights2[0][0][0].GetMask();
	class Line *msk1line = m_Weights2[0][1][0].GetMask();
	class Line *msk2line = m_Weights2[1][0][0].GetMask();
	class Line *msk3line = m_Weights2[1][1][0].GetMask();
	
	for(i = 0;i < w;i++) {
	  mskdest->At(i) = msk1line->Get(i) * 0.5f + msk2line->Get(i) * 0.5f + msk3line->Get(i) * 0.25f;
	}
      }
      //
      for(y = 0;y < 8;y++) {
	for(x = 0;x < 8;x++) {
	  // After warmup, all masks guaranteed valid (never nullptr)
	  class Line *refmask = m_Weights1[x][y][c].GetMask();
	  class Line *refline = m_Weights1[x][y][c].GetCoeff();
	  class Line *dstmask = m_Weights2[x][y][c].GetMask();
	  class Line *dstline = m_Weights2[x][y][c].GetCoeff();
	  
	  if (logme > 2) {
	    m_MaskImgs1[x][y][c].WriteLine(refmask);
	    m_MaskImgs2[x][y][c].WriteLine(dstmask);
	  }
	  
	  // Process this band (no nullptr checks needed after warmup)
	  {
#ifdef EXTENDED_FILTER
	    if (c == 0 && x == 0 && y == 0) {
	      // The LL band of the Y-channel
	      MeasureInBand(m_Weights1[0][0][0].GetLowpass(0,0),m_Weights2[0][0][0].GetLowpass(0,0),
			    refmask,dstmask,w,logme,luma_tbl[0]/90.0f,&errorline);
	      MeasureInBand(m_Weights1[0][0][0].GetLowpass(1,0),m_Weights2[0][0][0].GetLowpass(1,0),
			    refmask,dstmask,w,logme,luma_tbl[0]/15.0f,&errorline);
	      MeasureInBand(m_Weights1[0][0][0].GetLowpass(0,1),m_Weights2[0][0][0].GetLowpass(0,1),
			    refmask,dstmask,w,logme,luma_tbl[0]/15.0f,&errorline);
	      MeasureInBand(m_Weights1[0][0][0].GetLowpass(1,1),m_Weights2[0][0][0].GetLowpass(1,1),
			    refmask,dstmask,w,logme,luma_tbl[0]/20.0f,&errorline);
	    } else 
#endif
	      MeasureInBand(refline,dstline,refmask,dstmask,w,logme,1.0,&errorline);
	    collected = true;
	  }
	  if (collected && logme > 2) {
	    assert(refline);
	    m_DCTDiff[x][y][c].WriteLine(refline);
#ifdef EXTENDED_FILTER
	    if (x == 0 && y == 0 && c == 0) {
	      for(int j = 0; j < 2; j++) {
		for(int ii = 0; ii < 2; ii++) {
		  m_LowPasses1[ii][j].WriteLine(m_Weights1[0][0][0].GetLowpass(ii,j));
		  m_LowPasses2[ii][j].WriteLine(m_Weights2[0][0][0].GetLowpass(ii,j));
		}
	      }
	    }
#endif	
	  }
	}
      }
    }
    //
    // One line done. If errors got collected, update the error map.
    if (collected) {
      class Line *refsaliency = nullptr;
      class Line *dstsaliency = nullptr;
      //
      if (sref && sdist) {
	refsaliency = sref->ReadNextPFMLine();
	dstsaliency = sdist->ReadNextPFMLine();
      }
      //
      for(i = 0;i < w;i++) {
#if defined(WEIGHT_DELTA_E)
	float v = 1.0 - errorline.Get(i); 
	float t = 0;
	float a,c1,c2;
	float a2,c12,c22;
	//
	a  = m_Weights1[0][0][0].GetOriginal()->Get(i);
	a2 = m_Weights2[0][0][0].GetOriginal()->Get(i);
	if (comps == 3) {
	  // Note that the order is reversed.
	  c1  = m_Weights1[0][0][2].GetOriginal()->Get(i);
	  c12 = m_Weights2[0][0][2].GetOriginal()->Get(i);
	  c2  = m_Weights1[0][0][1].GetOriginal()->Get(i);
	  c22 = m_Weights2[0][0][1].GetOriginal()->Get(i);
	} else {
	  // Pure grey-scale otherwise.
	  c1  = 0;
	  c12 = 0;
	  c2  = 0;
	  c22 = 0;
	}
	float deltaE = DeltaE(a,c1,c2,a2,c12,c22);
	assert(deltaE >= 0.0);
	t     = deltaE * v;
#elif defined(WEIGHT_MSE)
	float v = 1.0 - errorline.Get(i); 
	float t = 0;
	for(c = 0;c < comps;c++) {
	  float org = m_Weights1[0][0][c].GetOriginal()->Get(i);
	  float dst = m_Weights2[0][0][c].GetOriginal()->Get(i);
	  t         += 64.0 * (org - dst)*(org - dst) * v;
	}
#else
	float t = 1.0f - errorline.Get(i);
	assert(t >= 0.0f && t <= 1.0f);
#endif
	// Include the saliency if we have it.
	if (refsaliency && dstsaliency) {
	  float p1,p2;
	  // t is the probability of finding an error. An error is found if *either* the error is
	  // salient in the reference or the distorted image, which is means that the error is only
	  // not visible if it is neither salient in reference or distorted image.
	  p1 = refsaliency->At(i + (offset >> 1));
	  p2 = dstsaliency->At(i + (offset >> 1));
	  t *= (1.0f - (1.0f - p1) * (1.0f - p2));
	}
	error   += t;
	errorline.At(i) = t;
      }	
      l++;
      if (logme) {
	m_ErrorMap.WriteLine(&errorline);
	printf(".");
	fflush(stdout);
      }
    }
  }
  //
  if (logme)
    printf("\n");
  // Take the average, then the logarithm.
  return -20.0f * logf(error / (w * h * comps)) / logf(10.0f);
}
///

#endif // HWY_ONCE
