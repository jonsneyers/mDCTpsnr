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
** $Id: colortransformer.cpp,v 1.5 2025/07/08 14:35:22 thor Exp $
**
*/

/// Includes
#include "options.hpp"
#include "ctrafo/colortransformer.hpp"
#include "dct/line.hpp"
#include "global/exceptions.hpp"
#include <cmath>
#include <cassert>
#include <algorithm>  // for std::min, std::max
///

/// PSI_Lookup
# define PSI_RAW(y) ((y) / ((y) + 12.6f*powf((y), 0.63f )))
#ifdef KNEE_VALUE
# define KNEE_THRES 12
# define KNEE_SLOPE (PSI_RAW(KNEE_THRES)/KNEE_THRES)
# define PSI_LOOKUP(y) ((((y) <= KNEE_THRES)?((y)*KNEE_SLOPE):PSI_RAW(y)))
#else
# define PSI_LOOKUP(y) ((((y) <= 0.0f)?(0.0f):PSI_RAW(y)))
#endif
///

/// ColorTransformer::ColorTransformer
ColorTransformer::ColorTransformer(void)
  : m_pdSRGBLinear(nullptr), m_pdLMS(nullptr)
{
  
}
///

/// ColorTransformer::~ColorTransformer
ColorTransformer::~ColorTransformer(void)
{
  delete[] m_pdSRGBLinear;
  delete[] m_pdLMS;
}
///

/// ColorTransformer::CreateLMSLookup
// Create the LMS lookup table.
void ColorTransformer::CreateLMSLookup(uint32_t scale)
{
  if (m_pdLMS == nullptr) {
    m_pdLMS = new float[scale];
    for(uint32_t i = 0; i < scale; i++) {
      // Inline LMS transfer: pow(in, 0.43)
      m_pdLMS[i] = powf(float(i) / (scale - 1), 0.43f);
    }
  }
}
///

/// ColorTransformer::CreateSRGBLinearLookup
// Create the sRGB to linear RGB lookup table based on bit depth
void ColorTransformer::CreateSRGBLinearLookup(uint8_t bits)
{
  if (m_pdSRGBLinear == nullptr && bits > 0) {
    uint32_t scale = 1UL << bits;
    m_pdSRGBLinear = new float[scale];
    
    for (uint32_t i = 0; i < scale; i++) {
      float normalized = float(i) / (scale - 1);
      // sRGB to linear RGB conversion
      if (normalized <= 0.04045f) {
        m_pdSRGBLinear[i] = normalized / 12.92f;
      } else {
        m_pdSRGBLinear[i] = powf((normalized + 0.055f) / 1.055f, 2.4f);
      }
    }
  }
}
///

/// ColorTransformer::ForwardsTransform
// Run a forwards transformation of the data
// in case we have three or more components. 
// Run it only on the first three.
// Input: Integer values (0..2^bits-1) stored as floats
// Output: Linearized and color-transformed floats
void ColorTransformer::ForwardsTransform(class Line *red,class Line *green,class Line *blue,uint8_t bits)
{
  uint32_t xp,width;
#ifdef USE_PSI
  float pmax = PSI_LOOKUP(100.0f);
#endif
#if defined (KNEE_VALUE) && !defined(USE_PSI)
  float kneeslope =  (1.055f * powf( 0.0031308f, float(OUTPUT_GAMMA) * 1.0f/2.4f ) - 0.055f) / 0.0031308f;
#endif
  //
  // Create sRGB linearization lookup table based on bit depth
  CreateSRGBLinearLookup(bits);
  //
#if defined(ITP)
  // Precision of 16 bit should be hopefully sufficient for LMS lookup.
  const uint32_t lmsscale = 1UL << 16;
  CreateLMSLookup(lmsscale);
#endif
  //
  width  = red->LengthOf();
  //
  for(xp = 0;xp<width;xp++) {
    // Get color values as integers (stored as floats)
    uint32_t r_int = uint32_t(red->Get(xp));
    uint32_t g_int = uint32_t(green->Get(xp));
    uint32_t b_int = uint32_t(blue->Get(xp));
    
    // Linearize via lookup table (sRGB to linear RGB)
    float r = m_pdSRGBLinear[r_int];
    float g = m_pdSRGBLinear[g_int];
    float b = m_pdSRGBLinear[b_int];
    //
#if defined(OPPONENT_COLOR)
    //
    float x = 0.4124564 * r + 0.2126729 * g + 0.0193339 * b;
    float y = 0.3575761 * r + 0.7151522 * g + 0.1191920 * b;
    float z = 0.1804375 * r + 0.0721750 * g + 0.9503041 * b;
    //
    // Thrid, convert from XYZ to ac1c2 (NOTE: The more precise figures come from the ICAM software)
    float a  =  0.2787336 * x + 0.7218031 * y - 0.1065520 * z;
    float c1 = -0.4487736 * x + 0.2898056 * y - 0.0771569 * z;
    float c2 =  0.0859513 * x - 0.5899859 * y + 0.5011089 * z;
    //
#ifdef KNEE_VALUE
    a = (a <= 0.0031308f ? a * kneeslope : 1.055f * powf(a, OUTPUT_GAMMA * 1.f/2.4f ) - 0.055f);
    /*
    cb = (cb <= 0.0031308f ? cb * 12.92f : 1.055f * powf( cb, OUTPUT_GAMMA * 1.f/2.4f ) - 0.055f);
    cr = (cr <= 0.0031308f ? cr * 12.92f : 1.055f * powf( cr, OUTPUT_GAMMA * 1.f/2.4f ) - 0.055f);
    */
#else
    // Simplified without knee - use powf for float precision
    a = powf(a,OUTPUT_GAMMA / 2.4f);
    /*
    cb = powf(cb,OUTPUT_GAMMA / 2.4f);
    cr = powf(cr,OUTPUT_GAMMA / 2.4f);
    */
#endif
    // Put the stuff back.
    red->Put(xp,a);
    green->Put(xp,c2); // c2 is yellow-blow, similar to cb
    blue->Put(xp,c1);  // c1 is red-green, similar to cr
#elif defined(ITP)
    float x,y,z;
    float l,m,s,i,t,p;
    // r, g, b are already linearized via m_pdSRGBLinear lookup table above
    //
    // Note that we need the D65 illuminant, which is the
    // illuminant of sRGB / Rec. 709.
    x  = 0.412453 * r + 0.357580 * g + 0.180423 * b;
    y  = 0.212671 * r + 0.715160 * g + 0.072169 * b;
    z  = 0.019334 * r + 0.119193 * g + 0.950227 * b;
    //
    // Convert from xyz to LMS. Yes, this could be done in one
    // step, but for simplicity...
    l  =  0.4002 * x + 0.7075 * y - 0.0807 * z;
    m  = -0.2280 * x + 1.1500 * y + 0.0612 * z;
    s  = 0.9184  * z;
    //
    // Nonlinear transformation into primed coordinates.
    l *= lmsscale;
    m *= lmsscale;
    s *= lmsscale;
    l  = (l >= 0.0f)?(m_pdLMS[uint32_t(l)]):(-m_pdLMS[uint32_t(-l)]);
    m  = (m >= 0.0f)?(m_pdLMS[uint32_t(m)]):(-m_pdLMS[uint32_t(-m)]);
    s  = (s >= 0.0f)?(m_pdLMS[uint32_t(s)]):(-m_pdLMS[uint32_t(-s)]);
    //
    // Transfer again from lms to IPT
    i  = 0.4    * l + 0.4    * m + 0.2    * s;
    p  = 4.4550 * l - 4.8510 * m + 0.3960 * s;
    t  = 0.8056 * l + 0.3572 * m - 1.1628 * s;
    //
    // Put the stuff back.
    red->Put(xp,i);
    green->Put(xp,p);
    blue->Put(xp,t);
#elif defined(LUV)
    //
    // First convert sRGB->XYZ
    float x,y,z;
    float xn,yn,zn,n;
    float l,u,v;
    float un,vn;
    //
    // r, g, b are already linearized via m_pdSRGBLinear lookup table above
    //
    // Note that we need the D65 illuminant, which is the
    // illuminant of sRGB / Rec. 709.
    x  = 0.412453 * r + 0.357580 * g + 0.180423 * b;
    y  = 0.212671 * r + 0.715160 * g + 0.072169 * b;
    z  = 0.019334 * r + 0.119193 * g + 0.950227 * b;
    //
    // Initialize with the D65 white-point.
    xn = 0.95056f;
    yn = 1.0f;
    zn = 1.089050f;
    //
    // Compute L*,u*,v*
    l  = 116.0f * powf(y, 1.0f/3.0f) - 16.0f;
    n  = (x + 15.0f*y + 3.0f * z);
    // Color coordinates for black are arbitrary. Set them zero.
    u  = (n != 0.0f)?(4.0f * x / n):(0.0f);
    v  = (n != 0.0f)?(9.0f * y / n):(0.0f);
    n  = (xn + 15.0f * yn + 3.0f * zn);
    un = 4.0f * xn / n;
    vn = 4.0f * yn / n;
    u  = 13.0f * l * (u - un);
    v  = 13.0f * l * (v - vn);
    
    assert(!isnan(l) && !isnan(u) && !isnan(v));
    
    //
    // Output scaling.
    //
    red->Put(xp  ,l / 100.0f);
    green->Put(xp,u / 100.0f);
    blue->Put(xp ,v / 100.0f);	
#elif defined(YOZ)
    // This is the Ahumada & Watson YOZ color space.
    // First convert sRGB->XYZ
    float x,y,z,o;
    // r, g, b are already linearized via m_pdSRGBLinear lookup table above
    //
    // Note that we need the D65 illuminant, which is the
    // illuminant of sRGB / Rec. 709.
    x  = 0.412424f * r + 0.212656f * g + 0.0193324f * b;
    y  = 0.357579f * r + 0.715158f * g + 0.119193f  * b;
    z  = 0.180464f * r + 0.0721856f* g + 0.950444f  * b;
    o  = 0.47f * x - 0.37f * y - 0.10f * z;
    //
    red->Put(xp,y);
    green->Put(xp,o);
    blue->Put(xp,z);
    //
#elif defined(LINEAR)
    float yv,cb,cr;
    // r, g, b are already linearized via lookup table above
    //
    // Transform to linear YCbCr, note that Y is in the *linear* domain.
    yv =  0.299f*r   + 0.587f  *g + 0.114f*b;
    cb = -0.16875f*r  - 0.33126f*g  + 0.5f*b;
    cr =  0.5f*r      - 0.41869f*g  - 0.08131f*b;
    //
    // Transform back to the nonlinear space.
#ifdef USE_PSI
    yv = PSI_LOOKUP(100.0f*yv) / pmax;
    /*
    cb = PSI_LOOKUP(100.0*cb) / pmax;
    cr = PSI_LOOKUP(100.0*cr) / pmax;
    */
#else
#ifdef KNEE_VALUE
    yv = (yv <= 0.0031308f ? yv * kneeslope : 1.055f * powf( yv, float(OUTPUT_GAMMA) * 1.0f/2.4f ) - 0.055f);
    /*
    cb = (cb <= 0.0031308f ? cb * 12.92f : 1.055f * powf( cb, OUTPUT_GAMMA * 1.f/2.4f ) - 0.055f);
    cr = (cr <= 0.0031308f ? cr * 12.92f : 1.055f * powf( cr, OUTPUT_GAMMA * 1.f/2.4f ) - 0.055f);
    */
#else
    // Simplified without knee - use powf for float precision
    yv = powf(yv,OUTPUT_GAMMA / 2.4f);
    /*
    cb = powf(cb,OUTPUT_GAMMA / 2.4f);
    cr = powf(cr,OUTPUT_GAMMA / 2.4f);
    */
#endif
#endif
    // Put the stuff back.
    red->Put(xp,yv);
    green->Put(xp,cb);
    blue->Put(xp,cr);
#else
    //
    // Regular YCbCr using Rec.601, also as used by JPEG-1.
    float yv,cb,cr;
    // Run the RGB->YCbCr transformation
    yv =  0.299f*r   + 0.587f  *g + 0.114f*b;
    cb = -0.16875f*r - 0.33126f*g + 0.5f*b;
    cr =  0.5f*r     - 0.41869f*g - 0.08131f*b;
    // Put the stuff back.
    red->Put(xp,yv);
    green->Put(xp,cb);
    blue->Put(xp,cr);
#endif
  }
}
///
