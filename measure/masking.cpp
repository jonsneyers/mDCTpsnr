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
#include "std/math.hpp"
#include "std/string.hpp"
#include "dct/line.hpp"
///

/// Masking::Masking
Masking::Masking(void)
  : m_pMask(NULL), m_pOutput(NULL), m_pOriginal(NULL), 
    m_pMapped(NULL), m_ulY(0), m_dMaskingSlope(1.0), m_dVisibility(1.0),
    m_bPostFilter(false)
  // Slope and masking exponent will be installed over.
{
  int x,y;
  
  memset(m_pInput   ,0,sizeof(m_pInput));
  memset(m_pBuffer  ,0,sizeof(m_pBuffer));
  memset(m_pAdded   ,0,sizeof(m_pAdded));
  memset(m_pFiltered,0,sizeof(m_pFiltered));

  // Initalize the window function. This is a Hamming window here.
  double xz = double(MaskSize - 1) * 0.5;
  double yz = double(MaskSize - 1) * 0.5;
  double xf = double(MaskSize + 1) * 0.5;
  double yf = double(MaskSize + 1) * 0.5;
  double norm = 0.0;
  double nsqr;

  for(y = 0;y < MaskSize;y++) {
    double gy = cos((y - yz) / yf * M_PI * 0.5);
    for(x = 0;x < MaskSize;x++) {
      double g = cos((x - xz) / xf * M_PI * 0.5) * gy;
      norm     += g;
      m_Window[x][y] = g;
    }
    m_Mask[y] = gy;
  }

  norm /= MaskSize * MaskSize;
  nsqr  = sqrt(norm);

  for(y = 0;y < MaskSize;y++) {
    m_Mask[y] /= nsqr;
    for(x = 0;x < MaskSize;x++) {
      m_Window[x][y] /= norm;
    }
  } 
}
///

/// Masking::~Masking
Masking::~Masking(void)
{
  int i,j;

  delete m_pMask;
  for(i = 0;i < MaskSize;i++) {
    delete m_pInput[i];
  }
  delete m_pMapped;
  delete m_pOriginal;
  delete m_pOutput;
  for(i = 0;i < MaskSize;i++) {
    delete m_pAdded[i];
    delete m_pBuffer[i];
  }
  for(j = 0;j < 2;j++) {
    for(i = 0;i < 2;i++) {
      delete m_pFiltered[i][j];
    }
  }
}
///

/// Masking::EnablePostFilter
// Install the post-filter
void Masking::EnablePostFilter(void)
{
  int i;
  m_bPostFilter = true;
  m_NormLo      = 0.0;
  m_NormHi      = 0.0;

  for(i=0;i < MaskSize;i++) {
    double xp = (i - ((MaskSize - 1) >> 1)) / double((MaskSize - 1) >> 1);
    double lo = 1.0;
    double hi = sin(xp * M_PI * 0.5);
    m_LowFilter[i] = lo;
    m_HiFilter[i]  = hi;
    m_NormHi      += hi * hi;
  }
  m_NormLo = 1.0 / MaskSize;
  m_NormHi = 1.0 / sqrt(m_NormHi);
}
///

/// Masking::PushOriginal
// Push an original line into the line buffer
void Masking::PushOriginal(class Line *line)
{
  if (m_pBuffer[m_ulY] == NULL) {
    m_pBuffer[m_ulY] = new class Line(line->LengthOf());
  }
  m_pBuffer[m_ulY]->Swap(*line);
}
///

/// Masking::PushLine
// Provide a new input line to the masking computation.
void Masking::PushLine(class Line *line)
{
  const DOUBLE vis     = m_dVisibility;
#ifdef DALY_MASKING
  const DOUBLE k1      = 0.01528670024;
  const DOUBLE k2      = 392.4980478;
  DOUBLE s  = k1 * pow(k2,m_dMaskingSlope);
  s         = 1.0 / sqrt(sqrt(s * s * s * s + 1.0));
#endif
  int i,w   = line->LengthOf();
  int width = w - MaskSize;
  //
  assert(width > 0);
  //
  // First check whether we have a corresponding input line available. If
  // not, create one now.
  if (m_pInput[m_ulY] == NULL) {
    m_pInput[m_ulY] = new class Line(w);
  }
  if (m_pMapped == NULL) {
    m_pMapped = new class Line(w);
  }
  if (m_pAdded[m_ulY] == NULL) {
    m_pAdded[m_ulY] = new class Line(w);
  }
  //
  // Swap the contents of the input line into our buffer.
  m_pInput[m_ulY]->Swap(*line);
  //
  // Convert to nonlinear space for the l^p norm - precompute to
  // avoid double computation.
  for(i = 0;i < w;i++) {
    //m_pMapped[m_ulY]->At(i) = sqrt(fabs(m_pInput[m_ulY]->Get(i)) / nominal);
    //
    // This follows closely the masking formula found in the VDP sources.
    DOUBLE v = fabs(m_pInput[m_ulY]->Get(i)) * vis;
#ifdef DALY_MASKING
    v = k1 * pow(k2 * v,m_dMaskingSlope);
    v = v * v * v * v;
    v = s * sqrt(sqrt(v + 1.0));
#else
    v = pow(v,m_dMaskingSlope);
#endif
    m_pMapped->At(i) = v; 
  }
  //
  // Already perform the column sum.
  for(i = 0;i < width;i++) {
    DOUBLE sum = 0.0;
    int x;
    for(x = 0;x < MaskSize;x++) {
      sum += m_Mask[x] * m_pMapped->Get(i+x);
    }
    m_pAdded[m_ulY]->At(i) = sum;
  }
  //
  // If the input buffer is full, compute the new mask.
  if (m_ulY == MaskSize-1) {
    class Line *top,*add,*buf;
    int y;
    //
    ComputeMask();
    //
#ifdef EXTENDED_FILTER
    if (m_bPostFilter) {
      ComputeLowpass();
    }
#endif
    //
    // Also move the proper buffer line to the output.
    if (m_pBuffer[MaskSize >> 1]) {
      if (m_pOriginal == NULL)
	m_pOriginal = new Line(width);
      for(i = 0;i < width;i++)
	m_pOriginal->At(i) = m_pBuffer[MaskSize >> 1]->At(i + (MaskSize >> 1));
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
// two.
void Masking::ComputeLowpass(void)
{
  double sll  = m_dVisibility * m_NormLo * m_NormLo;
  double slh  = m_dVisibility * m_NormLo * m_NormHi;
  double shh  = m_dVisibility * m_NormHi * m_NormHi;
  int width   = m_pInput[0]->LengthOf() - MaskSize;
  int i,j,k;

  for(j = 0;j < 2;j++) {
    for(i = 0;i < 2;i++) {
      if (m_pFiltered[i][j] == NULL)
	m_pFiltered[i][j] = new Line(width);
      m_pFiltered[i][j]->Zero();
    }
  }

  for(k = 0;k < width;k++) {
    for(j = 0;j < MaskSize;j++) {
      for(i = 0;i < MaskSize;i++) {
	m_pFiltered[0][0]->At(k) += m_pInput[i]->Get(j+k) * m_LowFilter[i] * m_LowFilter[j] * sll;
	m_pFiltered[1][0]->At(k) += m_pInput[i]->Get(j+k) * m_LowFilter[i] * m_HiFilter[j]  * slh;
	m_pFiltered[0][1]->At(k) += m_pInput[i]->Get(j+k) * m_HiFilter[i]  * m_LowFilter[j] * slh;
	m_pFiltered[1][1]->At(k) += m_pInput[i]->Get(j+k) * m_HiFilter[i]  * m_HiFilter[j]  * shh;
      }
    }
  }
}
///

/// Masking::ComputeMask
// Compute the masking strength for the buffered lines, create a new buffered
// output line.
void Masking::ComputeMask(void)
{
  int i,y;
  int width = m_pInput[0]->LengthOf() - MaskSize;
  DOUBLE sum;
  const DOUBLE base    = BASE_VISIBILITY;  // base visibility.
  //
  //
  assert(width > 0);
  if (m_pMask == NULL) {
    m_pMask = new Line(width);
  }
  if (m_pOutput == NULL) {
    m_pOutput = new Line(width);
  }
  // For the computation of the masking strength, follow the idea of Taubman and
  // compute a (weighted) l^p norm of the coefficients surrounding the target
  // coefficient. p is simply 0.5.
  // NOTE: This is the slow one. Any improvement here improves the total
  // running speed big time!
  for(i = 0;i < width;i++) {
    sum = 0.0;
    for(y = 0;y < MaskSize;y++) {
      sum += m_Mask[y] * m_pAdded[y]->Get(i);
    }
#if NO_BASE_VISIBILITY
    m_pMask->At(i)   = 1.0 / (sum / (MaskSize * MaskSize));
#else
    m_pMask->At(i)   = base / (base + sum / (MaskSize * MaskSize));
#endif
    // Also copy to the output using the same offset.
    m_pOutput->At(i) = m_pInput[MaskSize >> 1]->Get(i + (MaskSize >> 1)) * m_dVisibility;
  }
}
///

