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
#include <cstdint>
#include <cstdio>
///

/// Forwards
class Line;
///

/// Class ImgWriter
// This simple class writes out images given a base name and additional parameters.
class ImgWriter {
  //
  uint32_t  m_ulWidth;
  uint32_t  m_ulHeight;
  uint32_t  m_ulY;
  //
  // First scale by scale, then add offset.
  float m_dOffset;
  float m_dScale;
  //
  FILE *m_pFile;
  //
  char *m_pcName;
  //
public:
  ImgWriter(void);
  ~ImgWriter(void);
  //
  void OpenPGM(uint32_t width,uint32_t height,float scale,float offset,const char *basename,...);
  //
  // Write out a single line.
  void WriteLine(const class Line *line);
  //
};
///

///
#endif
