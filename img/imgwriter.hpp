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

#ifndef IMG_IMGWRITER_HPP
#define IMG_IMGWRITER_HPP

/// Includes
#include "global/types.hpp"
#include "std/stdio.hpp"
///

/// Forwards
class Line;
///

/// Class ImgWriter
// This simple class writes out images given a base name and additional parameters.
class ImgWriter {
  //
  ULONG  m_ulWidth;
  ULONG  m_ulHeight;
  ULONG  m_ulY;
  //
  // First scale by scale, then add offset.
  DOUBLE m_dOffset;
  DOUBLE m_dScale;
  //
  FILE *m_pFile;
  //
  char *m_pcName;
  //
public:
  ImgWriter(void);
  ~ImgWriter(void);
  //
  void OpenPGM(ULONG width,ULONG height,double scale,double offset,const char *basename,...);
  //
  // Write out a single line.
  void WriteLine(const class Line *line);
  //
};
///

///
#endif
