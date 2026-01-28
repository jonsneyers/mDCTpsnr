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
#include "std/stdio.hpp"
#include "std/math.hpp"
#include "measure/ediff.hpp"
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
static const DOUBLE dct_scale[] = {
  1.0, 1.387039845, 1.306562965, 1.175875602,
  1.0, 0.785694958, 0.541196100, 0.275899379
};

//
// Additional correction factors that make the transformation
// equal-frequency dependent. It isn't currently due to the window
// function. What we see here is the l^2 norm squared under the
// transformation on an input that is given by a constant-frequency
// DCT pattern.
static const DOUBLE norms[8][8] = {
  {1,0.742298,1.22273,1.24434,1.2528,1.26126,1.28286,1.76329},
  {0.742298,0.439822,0.724487,0.737285,0.742298,0.747311,0.760109,1.04477},
  {1.22273,0.724487,1.19339,1.21448,1.22273,1.23099,1.25207,1.72098},
  {1.24434,0.737285,1.21448,1.23593,1.24434,1.25274,1.27419,1.75139},
  {1.2528,0.742298,1.22273,1.24434,1.2528,1.26126,1.28286,1.76329},
  {1.26126,0.747311,1.23099,1.25274,1.26126,1.26977,1.29152,1.7752},
  {1.28286,0.760109,1.25207,1.27419,1.28286,1.29152,1.31364,1.8056},
  {1.76329,1.04477,1.72098,1.75139,1.76329,1.7752,1.8056,2.48181},
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
	    class Line *refline;
	    class Line *dstline;
	    //
	    if ((refline = m_pReference->GetDCTBand(x,y,c))) {
	      if (logme > 1) m_DCTImgs1[x][y][c].WriteLine(refline);
	      m_Weights1[x][y][c].PushLine(refline);
	    }
	    if ((dstline = m_pDistorted->GetDCTBand(x,y,c))) {
	    if (logme > 1) m_DCTImgs2[x][y][c].WriteLine(dstline);
	    m_Weights2[x][y][c].PushLine(dstline);
	    }
	  }
	  cnt++;
	}
      }
    }
  }
}
///

/// Pooling::MeasureInBand
// Measure the visibility in a band
// and adjust the error map.
void Pooling::MeasureInBand(class Line *refline,class Line *dstline,
			    const class Line *refmask,const class Line *dstmask,
			    int w,int logme,DOUBLE visbase,
			    class Line *errorline)
{ 
  int i;

  for(i = 0;i < w;i++) {
    DOUBLE vis    = 1.0;
    DOUBLE visref = 1.0;
    DOUBLE visdst = 1.0;
    DOUBLE err;
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
    err = fabs(err); // DCTPSNR uses here err*err, but the difference gets visible, not its square
    err = err * vis;
    assert(err >= 0.0);
    if (err > 0.0) {
      // The visibility threshold is normalized such that
      // 1.0 is approximately 50% detection probabililty.
      // However, 1.0 is the "noticability" threshold, thus
      // upscale.
#ifdef AHUMADA
      err = exp(AHUMADA_FACTOR * -pow(DETECT_THREAS * err,AHUMADA_EXPONENT)); // 1.0 - detection probability.
#else
      err = exp(-pow(2.0 * err,3.5)); // 1.0 - detection probability.
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
DOUBLE Pooling::Measure(class Image *ref ,class Image *dist,
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
  DOUBLE error = 0.0;
  //
  //
  m_iComponents = comps;
  m_iLogging    = logme;
  m_pReference  = ref;
  m_pDistorted  = dist;
  //
  // If saliency maps are defined, read the first lines to correct for the offsets.
  if (sref && sdist) {
    int i;
    for(i = 0;i < (offset >> 1);i++) {
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
	  DOUBLE offset,scale;
	  if (x == 0 && y == 0) {
	    offset = 0.0; 
	    scale  = 255.0/8.0; // show the average luminance, not the nominal low-pass
	  } else {
	    offset = 128.0;
	    scale  = 128.0;
	  }
	  // Correct for the wrong output scaling of the exponent.
	  scale /= (8.0 * dct_scale[x] * dct_scale[y] * sqrt(norms[x][y]));
	  m_DCTImgs1[x][y][c].OpenPGM(ref->WidthOf() - 8,ref->HeightOf() - 8,scale,offset,"dct1_%d_%d_%d.pgm",x,y,c);
	  m_DCTImgs2[x][y][c].OpenPGM(dist->WidthOf() - 8,dist->HeightOf() - 8,scale,offset,"dct2_%d_%d_%d.pgm",x,y,c);
	  m_DCTDiff[x][y][c].OpenPGM(w,h,128.0,128.0,"dctdiff_%d_%d_%d.pgm",x,y,c);
	  m_MaskImgs1[x][y][c].OpenPGM(w,h,255.0,0.0,"mask1_%d_%d_%d.pgm",x,y,c);
	  m_MaskImgs2[x][y][c].OpenPGM(w,h,255.0,0.0,"mask2_%d_%d_%d.pgm",x,y,c);
#ifdef EXTENDED_FILTER
	  if (x == 0 && y == 0 && c == 0) {
	    // DCT scale correction is already done in the visibility computation
	    // of the mask, do not repeat.
	    scale = luma_tbl[0] / 255.0;
	    m_LowPasses1[0][0].OpenPGM(w,h,scale * 255.0 / 8.0,0.0   ,"lowpass1_%d_%d.pgm",0,0);
	    m_LowPasses2[0][0].OpenPGM(w,h,scale * 255.0 / 8.0,0.0   ,"lowpass2_%d_%d.pgm",0,0);
	    m_LowPasses1[0][1].OpenPGM(w,h,scale * 128.0      ,128.0 ,"lowpass1_%d_%d.pgm",0,1);
	    m_LowPasses2[0][1].OpenPGM(w,h,scale * 128.0      ,128.0 ,"lowpass2_%d_%d.pgm",0,1);
	    m_LowPasses1[1][0].OpenPGM(w,h,scale * 128.0      ,128.0 ,"lowpass1_%d_%d.pgm",1,0);
	    m_LowPasses2[1][0].OpenPGM(w,h,scale * 128.0      ,128.0 ,"lowpass2_%d_%d.pgm",1,0);
	    m_LowPasses1[1][1].OpenPGM(w,h,scale * 128.0      ,128.0 ,"lowpass1_%d_%d.pgm",1,1);
	    m_LowPasses2[1][1].OpenPGM(w,h,scale * 128.0      ,128.0 ,"lowpass2_%d_%d.pgm",1,1);
	  }
#endif
	}
      }
      {
	ULONG wo = ref->WidthOf()  - 8;
	ULONG ho = ref->HeightOf() - 8;
	DOUBLE offset = (c == 0)?(0):(128);
	m_SrcImgs1[c].OpenPGM(wo,ho,255.0,offset,"refimg_%d.pgm",c);
	m_SrcImgs2[c].OpenPGM(wo,ho,255.0,offset,"dstimg_%d.pgm",c);
      }
    }
  }

  //
  // Setup the masking and visibility weights for the bands.
  for(c = 0;c < comps;c++) {
    for(y = 0;y < 8;y++) {
      for(x = 0;x < 8;x++) {
	DOUBLE vis,expon;
	int xs = x; // << 1;
	int ys = y; // << 1;
	if (xs > 7) xs = 7;
	if (ys > 7) ys = 7;
	//
	// First integrate the scaling of the DCT. Input must be scaled by 8, then by the
	// norm correction factor.
	vis = 1.0 / (8.0 * dct_scale[x] * dct_scale[y] * sqrt(norms[x][y]));
	//
	// The idea of the next scaling is to bring the visibility just to
	// the "just noticable difference", which is by an 8bpp input given
	// by the quantization coefficients of the JPEG compression. All other
	// coefficients are then scaled to the visibility threshold of the DC
	// luminance.
	vis *= 255.0 / 16.0;
	//
	// Further, include the base visibility of the frequency band.
	switch(c) {
	case 0: // Y
	  vis *= 1.0   * 16.0 / luma_tbl  [x + (y << 3)];
	  break;
	case 1: // Cb
#ifdef AHUMADA
	  vis *= 0.088 * 45.0 / cb_tbl    [xs + (ys << 3)];
#else
	  vis *= 0.088 * 16.0 / chroma_tbl[x + (y << 3)];
#endif
	  break;
	case 2: // Cr
#ifdef AHUMADA
	  vis *= 0.278 * 21.0 / cr_tbl    [xs + (ys << 3)];
#else
	  vis *= 0.278 * 16.0 / chroma_tbl[x + (y << 3)];
#endif
	  break;
	}
	//
	// Get the masking exponent. This is also frequency dependent.
	switch(x+y) {
	case 0:
	  expon = 0.7;
	  break;
	case 1:
	case 2:
	  expon = 0.8;
	  break;
	case 3:
	  expon = 0.9;
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
  // Now for the main loop.
  while(!ref->ImageDone()) {
    bool collected = false;
    ref->ReadNextLine();
    dist->ReadNextLine(); 
    //
    //
    // Reset the detection probability for this line.
    for(i = 0;i < w;i++) {
      errorline.At(i) = 1.0;
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
	    mskdest->At(i) = msk1line->Get(i) * 0.5 + msk2line->Get(i) * 0.5 + msk3line->Get(i) * 0.25;
	  }
	}
      }
      //
      {
	class Line *mskdest  = m_Weights2[0][0][0].GetMask();
	class Line *msk1line = m_Weights2[0][1][0].GetMask();
	class Line *msk2line = m_Weights2[1][0][0].GetMask();
	class Line *msk3line = m_Weights2[1][1][0].GetMask();
	//
	if (mskdest) {
	  for(i = 0;i < w;i++) {
	    mskdest->At(i) = msk1line->Get(i) * 0.5 + msk2line->Get(i) * 0.5 + msk3line->Get(i) * 0.25;
	  }
	}
      }
      //
      for(y = 0;y < 8;y++) {
	for(x = 0;x < 8;x++) {
	  class Line *refline = NULL,*refmask = NULL;
	  class Line *dstline = NULL,*dstmask = NULL;
	  //
	  if ((refmask = m_Weights1[x][y][c].GetMask())) {
	    if (logme > 2) m_MaskImgs1[x][y][c].WriteLine(refmask);
	    refline = m_Weights1[x][y][c].GetCoeff();
	  }
	  if ((dstmask = m_Weights2[x][y][c].GetMask())) {
	    if (logme > 2) m_MaskImgs2[x][y][c].WriteLine(dstmask);
	    dstline = m_Weights2[x][y][c].GetCoeff();
	  }
	  //
	  // Now iterate over the line of samples if we could collect
	  // data.
	  if (refmask && dstmask) {
#ifdef EXTENDED_FILTER
	    if (c == 0 && x == 0 && y == 0) {
	      // The LL band of the Y-channel
	      MeasureInBand(m_Weights1[0][0][0].GetLowpass(0,0),m_Weights2[0][0][0].GetLowpass(0,0),
			    refmask,dstmask,w,logme,luma_tbl[0]/90.0,&errorline);
	      MeasureInBand(m_Weights1[0][0][0].GetLowpass(1,0),m_Weights2[0][0][0].GetLowpass(1,0),
			    refmask,dstmask,w,logme,luma_tbl[0]/15.0,&errorline);
	      MeasureInBand(m_Weights1[0][0][0].GetLowpass(0,1),m_Weights2[0][0][0].GetLowpass(0,1),
			    refmask,dstmask,w,logme,luma_tbl[0]/15.0,&errorline);
	      MeasureInBand(m_Weights1[0][0][0].GetLowpass(1,1),m_Weights2[0][0][0].GetLowpass(1,1),
			    refmask,dstmask,w,logme,luma_tbl[0]/20.0,&errorline);
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
	      int i,j;
	      for(j = 0;j < 2;j++) {
		for(i = 0;i < 2;i++) {
		  m_LowPasses1[i][j].WriteLine(m_Weights1[0][0][0].GetLowpass(i,j));
		  m_LowPasses2[i][j].WriteLine(m_Weights2[0][0][0].GetLowpass(i,j));
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
      class Line *refsaliency = NULL;
      class Line *dstsaliency = NULL;
      //
      if (sref && sdist) {
	refsaliency = sref->ReadNextPFMLine();
	dstsaliency = sdist->ReadNextPFMLine();
      }
      //
      for(i = 0;i < w;i++) {
#if defined(WEIGHT_DELTA_E)
	DOUBLE v = 1.0 - errorline.Get(i); 
	DOUBLE t = 0;
	DOUBLE a,c1,c2;
	DOUBLE a2,c12,c22;
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
	DOUBLE deltaE = DeltaE(a,c1,c2,a2,c12,c22);
	assert(deltaE >= 0.0);
	t     = deltaE * v;
#elif defined(WEIGHT_MSE)
	DOUBLE v = 1.0 - errorline.Get(i); 
	DOUBLE t = 0;
	for(c = 0;c < comps;c++) {
	  DOUBLE org = m_Weights1[0][0][c].GetOriginal()->Get(i);
	  DOUBLE dst = m_Weights2[0][0][c].GetOriginal()->Get(i);
	  t         += 64.0 * (org - dst)*(org - dst) * v;
	}
#else
	DOUBLE t = 1.0 - errorline.Get(i);
	assert(t >= 0.0 && t <= 1.0);
#endif
	// Include the saliency if we have it.
	if (refsaliency && dstsaliency) {
	  DOUBLE p1,p2;
	  // t is the probability of finding an error. An error is found if *either* the error is
	  // salient in the reference or the distorted image, which is means that the error is only
	  // not visible if it is neither salient in reference or distorted image.
	  p1 = refsaliency->At(i + (offset >> 1));
	  p2 = dstsaliency->At(i + (offset >> 1));
	  t *= (1.0 - (1.0 - p1) * (1.0 - p2));
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
  return -20.0 * log(error / (w * h * comps)) / log(10.0);
}
///
