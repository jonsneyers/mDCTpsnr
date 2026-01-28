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
#include "std/math.hpp"
#include "std/assert.hpp"
///

/// PSI_Lookup
# define PSI_RAW(y) ((y) / ((y) + 12.6f*powf((y), 0.63f )))
#ifdef KNEE_VALUE
# define KNEE_THRES 12
# define KNEE_SLOPE (PSI_RAW(KNEE_THRES)/KNEE_THRES)
# define PSI_LOOKUP(y) ((((y) <= KNEE_THRES)?((y)*KNEE_SLOPE):PSI_RAW(y)))
#else
# define PSI_LOOKUP(y) ((((y) <= 0.0)?(0.0):PSI_RAW(y)))
#endif
///

/// ColorTransformer::ColorTransformer
ColorTransformer::ColorTransformer(void)
  : m_pdLookup(NULL), m_pdLMS(NULL)
{
  
}
///

/// ColorTransformer::~ColorTransformer
ColorTransformer::~ColorTransformer(void)
{
  delete[] m_pdLookup;
}
///

/// ColorTransformer::CreateLookup
// Create the lookup table.
void ColorTransformer::CreateLookup(ULONG scale)
{
  ULONG i;

  if (m_pdLookup == NULL) {
    m_pdLookup = new DOUBLE[scale];
    
    for(i = 0; i < scale;i++) {
      m_pdLookup[i] = sRGBTransfer(i,scale - 1);
    }
  }
}
///

/// ColorTransformer::CreateLMSLookup
// Create the LMS lookup table.
void ColorTransformer::CreateLMSLookup(ULONG scale)
{
  ULONG i;

  if (m_pdLMS == NULL) {
    m_pdLMS = new DOUBLE[scale];
    for(i = 0;i < scale;i++) {
      m_pdLMS[i] = LMSTransfer(DOUBLE(i) / (scale - 1));
    }
  }
}
///

/// ColorTransformer::ForwardsTransform
// Run a forwards transformation of the data
// in case we have three or more components. 
// Run it only on the first three.
void ColorTransformer::ForwardsTransform(class Line *red,class Line *green,class Line *blue,UBYTE)
{
  ULONG xp,width;
#ifdef USE_PSI
  double pmax = PSI_LOOKUP(100.0);
#endif
#if defined (KNEE_VALUE) && !defined(USE_PSI)
  double kneeslope =  (1.055 * pow( 0.0031308, OUTPUT_GAMMA * 1./2.4 ) - 0.055) / 0.0031308;
#endif
  //
  //
#if defined(LUV) || defined(ITP) || defined(YOZ)
  // Create the input lookup table to convert RGB to R'G'B'.
  const DOUBLE max     = 1.0;
  const ULONG scale    = 1UL << 12;
  const ULONG lmsscale = 1UL << 16;
  CreateLookup(scale);
  // Precision of 16 bit should be hopefully sufficient.
  CreateLMSLookup(lmsscale);
#endif
  //
  width  = red->LengthOf();
  //
  for(xp = 0;xp<width;xp++) {
    DOUBLE r,g,b;
    // Get color values. (The slow way).
    r  = red->Get(xp);
    g  = green->Get(xp);
    b  = blue->Get(xp);
    //
#if defined(OPPONENT_COLOR)
    // Linearized RGB
    r  = (r <= 0.04045 ? r / 12.92 : pow( (r + 0.055) / 1.055, 2.4));
    g  = (g <= 0.04045 ? g / 12.92 : pow( (g + 0.055) / 1.055, 2.4));
    b  = (b <= 0.04045 ? b / 12.92 : pow( (b + 0.055) / 1.055, 2.4));
    //
    double x = 0.4124564 * r + 0.2126729 * g + 0.0193339 * b;
    double y = 0.3575761 * r + 0.7151522 * g + 0.1191920 * b;
    double z = 0.1804375 * r + 0.0721750 * g + 0.9503041 * b;
    //
    // Thrid, convert from XYZ to ac1c2 (NOTE: The more precise figures come from the ICAM software)
    double a  =  0.2787336 * x + 0.7218031 * y - 0.1065520 * z;
    double c1 = -0.4487736 * x + 0.2898056 * y - 0.0771569 * z;
    double c2 =  0.0859513 * x - 0.5899859 * y + 0.5011089 * z;
    //
#ifdef KNEE_VALUE
    a = (a <= 0.0031308 ? a * kneeslope : 1.055 * pow(a, OUTPUT_GAMMA * 1./2.4 ) - 0.055);
    /*
    cb = (cb <= 0.0031308 ? cb * 12.92 : 1.055 * pow( cb, OUTPUT_GAMMA * 1./2.4 ) - 0.055);
    cr = (cr <= 0.0031308 ? cr * 12.92 : 1.055 * pow( cr, OUTPUT_GAMMA * 1./2.4 ) - 0.055);
    */
#else
    // Simplified without knee
    a = pow(a,OUTPUT_GAMMA / 2.4);
    /*
    cb = pow(cb,OUTPUT_GAMMA / 2.4);
    cr = pow(cr,OUTPUT_GAMMA / 2.4);
    */
#endif
    // Put the stuff back.
    red->Put(xp,a);
    green->Put(xp,c2); // c2 is yellow-blow, similar to cb
    blue->Put(xp,c1);  // c1 is red-green, similar to cr
#elif defined(ITP)
    DOUBLE x,y,z;
    DOUBLE l,m,s,i,t,p;
    //
    // Clip overflows.
    // First step: Transform sRGB to R'G'B'
    if (r <= 0.0) {
      r = 0;
    } else if (r >= max) {
      r = 1.0;
    } else {
      r = m_pdLookup[ULONG(r * scale)];
    }
    if (g <= 0.0) {
      g = 0;
    } else if (g >= max) {
      g = 1.0;
    } else {
      g = m_pdLookup[ULONG(g* scale)];
    }
    if (b <= 0.0) {
      b = 0;
    } else if (b >= max) {
      b = 1.0;
    } else {
      b = m_pdLookup[ULONG(b* scale)];
    }
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
    l  = (l >= 0.0)?(m_pdLMS[ULONG(l)]):(-m_pdLMS[ULONG(-l)]);
    m  = (m >= 0.0)?(m_pdLMS[ULONG(m)]):(-m_pdLMS[ULONG(-m)]);
    s  = (s >= 0.0)?(m_pdLMS[ULONG(s)]):(-m_pdLMS[ULONG(-s)]);
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
    DOUBLE x,y,z;
    DOUBLE xn,yn,zn,n;
    DOUBLE l,u,v;
    DOUBLE un,vn;
    //
    // Clip overflows.
    // First step: Transform sRGB to R'G'B'
    if (r <= 0.0) {
      r = 0;
    } else if (r >= max) {
      r = 1.0;
    } else {
      r = m_pdLookup[ULONG(r* scale)];
    }
    if (g <= 0.0) {
      g = 0;
    } else if (g >= max) {
      g = 1.0;
    } else {
      g = m_pdLookup[ULONG(g* scale)];
    }
    if (b <= 0.0) {
      b = 0;
    } else if (b >= max) {
      b = 1.0;
    } else {
      b = m_pdLookup[ULONG(b* scale)];
    }
    //
    // Note that we need the D65 illuminant, which is the
    // illuminant of sRGB / Rec. 709.
    x  = 0.412453 * r + 0.357580 * g + 0.180423 * b;
    y  = 0.212671 * r + 0.715160 * g + 0.072169 * b;
    z  = 0.019334 * r + 0.119193 * g + 0.950227 * b;
    //
    // Initialize with the D65 white-point.
    xn = 0.95056;
    yn = 1.0;
    zn = 1.089050;
    //
    // Compute L*,u*,v*
    l  = 116.0 * pow(y,1.0/3.0) - 16.0;
    n  = (x + 15*y + 3 * z);
    // Color coordinates for black are arbitrary. Set them zero.
    u  = (n != 0.0)?(4 * x / n):(0.0);
    v  = (n != 0.0)?(9 * y / n):(0.0);
    n  = (xn + 15 * yn + 3 * zn);
    un = 4 * xn / n;
    vn = 4 * yn / n;
    u  = 13 * l * (u - un);
    v  = 13 * l * (v - vn);
    
    assert(!isnan(l) && !isnan(u) && !isnan(v));
    
    //
    // Output scaling.
    //
    red->Put(xp  ,l / 100.0);
    green->Put(xp,u / 100.0);
    blue->Put(xp ,v / 100.0);	
#elif defined(YOZ)
    // This is the Ahumada & Watson YOZ color space.
    // First convert sRGB->XYZ
    DOUBLE x,y,z,o;
    //
    // Clip overflows.
    // First step: Transform sRGB to R'G'B'
    if (r <= 0.0) {
      r = 0;
    } else if (r >= max) {
      r = 1.0;
    } else {
      r = m_pdLookup[ULONG(r* scale)];
    }
    if (g <= 0.0) {
      g = 0;
    } else if (g >= max) {
      g = 1.0;
    } else {
      g = m_pdLookup[ULONG(g* scale)];
    }
    if (b <= 0.0) {
      b = 0;
    } else if (b >= max) {
      b = 1.0;
    } else {
      b = m_pdLookup[ULONG(b* scale)];
    }
    //
    // Note that we need the D65 illuminant, which is the
    // illuminant of sRGB / Rec. 709.
    x  = 0.412424 * r + 0.212656 * g + 0.0193324 * b;
    y  = 0.357579 * r + 0.715158 * g + 0.119193  * b;
    z  = 0.180464 * r + 0.0721856* g + 0.950444  * b;
    o  = 0.47 * x - 0.37 * y - 0.10 * z;
    //
    red->Put(xp,y);
    green->Put(xp,o);
    blue->Put(xp,z);
    //
#elif defined(LINEAR)
    DOUBLE yv,cb,cr;
    DOUBLE rp,gp,bp;
    //
    // Linearized RGB/YCbCr
    rp  = (r <= 0.04045 ? r / 12.92 : pow( (r + 0.055) / 1.055, 2.4));
    gp  = (g <= 0.04045 ? g / 12.92 : pow( (g + 0.055) / 1.055, 2.4));
    bp  = (b <= 0.04045 ? b / 12.92 : pow( (b + 0.055) / 1.055, 2.4));
    //
    // Transform to linear YCbCr, note that Y is in the *linear* domain.
    yv =  0.299*rp   + 0.587  *gp + 0.114*bp;
    cb = -0.16875*r  - 0.33126*g  + 0.5*b;
    cr =  0.5*r      - 0.41869*g  - 0.08131*b;
    //
    // Transform back to the nonlinear space.
#ifdef USE_PSI
    yv = PSI_LOOKUP(100.0*yv) / pmax;
    /*
    cb = PSI_LOOKUP(100.0*cb) / pmax;
    cr = PSI_LOOKUP(100.0*cr) / pmax;
    */
#else
#ifdef KNEE_VALUE
    yv = (yv <= 0.0031308 ? yv * kneeslope : 1.055 * pow( yv, OUTPUT_GAMMA * 1./2.4 ) - 0.055);
    /*
    cb = (cb <= 0.0031308 ? cb * 12.92 : 1.055 * pow( cb, OUTPUT_GAMMA * 1./2.4 ) - 0.055);
    cr = (cr <= 0.0031308 ? cr * 12.92 : 1.055 * pow( cr, OUTPUT_GAMMA * 1./2.4 ) - 0.055);
    */
#else
    // Simplified without knee
    yv = pow(yv,OUTPUT_GAMMA / 2.4);
    /*
    cb = pow(cb,OUTPUT_GAMMA / 2.4);
    cr = pow(cr,OUTPUT_GAMMA / 2.4);
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
    DOUBLE yv,cb,cr;
    // Run the RGB->YCbCr transformation
    yv =  0.299*r   + 0.587  *g + 0.114*b;
    cb = -0.16875*r - 0.33126*g + 0.5*b;
    cr =  0.5*r     - 0.41869*g - 0.08131*b;
    // Put the stuff back.
    red->Put(xp,yv);
    green->Put(xp,cb);
    blue->Put(xp,cr);
#endif
  }
}
///
