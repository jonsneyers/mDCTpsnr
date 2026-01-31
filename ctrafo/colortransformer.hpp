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
** $Id: colortransformer.hpp,v 1.4 2025/07/08 14:35:22 thor Exp $
**
*/

#ifndef CTRAFO_COLORTRANSFORMER_HPP
#define CTRAFO_COLORTRANSFORMER_HPP

/// Includes
#include <cmath>
#include <cstdint>
///

/// Forwards
class Line;
///

/// ColorTransformer
// This class describes a simplistic color transformer
// for RGB->YC_bC_r transformation.
class ColorTransformer {
  //
  // The lookup table for sRGB linearization (uint to float linear)
  float *m_pdSRGBLinear;
  //
  // The lookup table for the LMS transfer function.
  float *m_pdLMS;
  //
  // Create the LMS lookup table.
  void CreateLMSLookup(uint32_t scale);
  //
  // Create the sRGB linearization lookup table.
  void CreateSRGBLinearLookup(uint8_t bits);
  //
public:
  ColorTransformer(void);
  ~ColorTransformer(void);
  //
  // Forwards transform, i.e. RGB->YC_bC_r
  void ForwardsTransform(class Line *red,class Line *green, class Line *blue,uint8_t bitdepth);
};
///

///
#endif
