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
** $Id: masking.hpp,v 1.8 2025/07/08 14:35:22 thor Exp $
**
*/

#ifndef MEASURE_MASKING_HPP
#define MEASURE_MASKING_HPP

/// Includes
#include <cstdint>
///

/// Forwards
class Line;
///

/// class Masking
// This class computes the masking strength of a region of coefficients in the DCT domain.
class Masking {
  //
public:
  // Dimensions of the masking neighbourhood
  enum {
    MaskSize = 5
  };
  //
private:
  //
  // Line buffers (wrapper objects with minimal overhead)
  class Line    *m_pMask;
  class Line    *m_pOutput;
  class Line    *m_pOriginal;
  class Line    *m_pInput[MaskSize];
  class Line    *m_pBuffer[MaskSize];
  class Line    *m_pMapped;
  class Line    *m_pAdded[MaskSize];
  class Line    *m_pFiltered[2][2];
  //
  // Contiguous buffers for all data (no individual allocations)
  float        *m_pWideBuffer;    // For width+MaskSize lines
  float        *m_pNormalBuffer;  // For width lines
  uint32_t          m_ulWideStride;   // Stride for wide lines
  uint32_t          m_ulNormalStride; // Stride for normal lines
  uint32_t          m_ulWideWidth;    // Cached wide width (width + MaskSize)
  //
  // The one-dimensional window function (for masking)
  float         m_Mask[MaskSize];
  //
  // Current Y position.
  uint32_t          m_ulY;
  //
  // The masking exponent
  float         m_dMaskingSlope;
  //
  // The visibility of this band. Input is multiplied with this factor.
  float         m_dVisibility;
  //
  // Set in case another filter should be installed on the output.
  bool           m_bPostFilter;
  //
  // Filter coefficients for the post filter.
  float         m_LowFilter[MaskSize];
  float         m_HiFilter[MaskSize];
  //
  // Normalization constants for the low and high-passes
  float         m_NormLo,m_NormHi;
  //
  // Width of buffers (0 = not initialized yet)
  uint32_t          m_ulWidth;
  //
  // Cached filter coefficient products (precomputed for ComputeLowpass)
  float         m_FilterCoeffs_LL[MaskSize][MaskSize];
  float         m_FilterCoeffs_LH[MaskSize][MaskSize];
  float         m_FilterCoeffs_HL[MaskSize][MaskSize];
  float         m_FilterCoeffs_HH[MaskSize][MaskSize];
  bool           m_bCoeffsComputed;
  //
  // Initialize all buffers to the given width (called on first PushLine)
  void InitializeBuffers(uint32_t width) __attribute__((noinline));
  //
  // Compute the masking strength for the buffered lines, create a new buffered
  // output line.
  void ComputeMask(void);
  //
  // Compute the low and high-pass filter of the input, splitting one band into
  // two.
  void ComputeLowpass(void);
  //
public:
  Masking(void);
  ~Masking(void);
  //
  // Define the masking exponent for this band
  void SetMaskingExponent(float expon)
  {
    m_dMaskingSlope = expon;
  }
  //
  // Define the visibility of this band.
  void SetVisibility(float vis)
  {
    m_dVisibility  = vis;
  }
  //
  // Install the post-filter
  void EnablePostFilter(void);
  //
  // Provide a new input line to the masking computation.
  void PushLine(class Line *line);
  //
  // Push the original line into the buffer to make it
  // available on time together with the masked output.
  void PushOriginal(class Line *line);
  //
  // Return the current line describing the masking strenght, or nullptr if this is not
  // yet computed.
  class Line *GetMask(void) const
  {
    return m_pMask;
  }
  //
  // Return the line containing the pixels being masked. This output fits to the
  // line returned by GetMask.
  class Line *GetCoeff(void) const
  {
    return m_pOutput;
  }
  //
  // Return the original input line, offset'ed correctly to represent the same line
  // Returns nullptr - original line not needed in default mode
  class Line *GetOriginal(void) const
  {
    return nullptr;
  }
  //
  // Return the low-pass bands.
  class Line *GetLowpass(int x,int y) const
  {
    return m_pFiltered[x][y];
  }
};
///

///
#endif
