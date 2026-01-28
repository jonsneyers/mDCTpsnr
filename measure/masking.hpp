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
#include "global/types.hpp"
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
  // Line buffer for the lines containing the masking strength.
  class Line    *m_pMask;
  //
  // The output line: Shifted input, to be masked by the mask.
  class Line    *m_pOutput;
  //
  // The output line: Shifted input, unmasked.
  class Line    *m_pOriginal;
  //
  // Input buffer for the unmasked lines.
  class Line    *m_pInput[MaskSize];
  //
  // Input buffer for the original lines.
  class Line    *m_pBuffer[MaskSize];
  //
  // The scaled/nonlinear input lines.
  class Line    *m_pMapped;
  //
  // The line sum, multiplied with the window in X direction.
  class Line    *m_pAdded[MaskSize];
  //
  // The four bands of the low-pass filter output (LL,HL,LH and HH)
  class Line    *m_pFiltered[2][2];
  //
  // The window function for computing the mask.
  DOUBLE         m_Window[MaskSize][MaskSize];
  //
  // The one-dimensional window function.
  DOUBLE         m_Mask[MaskSize];
  //
  // Current Y position.
  ULONG          m_ulY;
  //
  // The masking exponent
  DOUBLE         m_dMaskingSlope;
  //
  // The visibility of this band. Input is multiplied with this factor.
  DOUBLE         m_dVisibility;
  //
  // Set in case another filter should be installed on the output.
  bool           m_bPostFilter;
  //
  // Filter coefficients for the post filter.
  DOUBLE         m_LowFilter[MaskSize];
  DOUBLE         m_HiFilter[MaskSize];
  //
  // Normalization constants for the low and high-passes
  DOUBLE         m_NormLo,m_NormHi;
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
  void SetMaskingExponent(DOUBLE expon)
  {
    m_dMaskingSlope = expon;
  }
  //
  // Define the visibility of this band.
  void SetVisibility(DOUBLE vis)
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
  // Return the current line describing the masking strenght, or NULL if this is not
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
  // as the masked line.
  class Line *GetOriginal(void) const
  {
    return m_pOriginal;
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
