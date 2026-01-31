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
#include "imgwriter.hpp"
#include <cstdio>
#include <cstdarg>
#include "dct/line.hpp"
#include <cerrno>
#include "global/exceptions.hpp"
///

/// ImgWriter::ImgWriter
ImgWriter::ImgWriter(void)
  : m_ulY(0), m_pFile(nullptr), m_pcName(nullptr)
{
}
///

/// ImgWriter::OpenPGM
void ImgWriter::OpenPGM(uint32_t w,uint32_t h,float scale,float offset,const char *basename,...)
{
  char buffer[256];
  va_list args;

  assert(m_pFile == nullptr);

  va_start(args,basename);
  vsnprintf(buffer,255,basename,args);
  va_end(args);

  m_ulWidth  = w;
  m_ulHeight = h;
  m_ulY      = 0;
  m_dScale   = scale;
  m_dOffset  = offset;
  m_pcName   = new char[strlen(buffer) + 1];
  strcpy(m_pcName,buffer);
}
///

/// ImgWriter::~ImgWriter
ImgWriter::~ImgWriter(void)
{
  if (m_pFile)
    fclose(m_pFile);
  delete[] m_pcName;
}
///

/// ImgWriter::WriteLine
void ImgWriter::WriteLine(const class Line *line)
{
  if (m_ulY == 0) { 
    m_pFile = fopen(m_pcName,"wb");
    if (m_pFile == nullptr) {
      ThrowIo("ImgWriter::ImgWriter","unable to open output file");
    }
    // Write header.
    fprintf(m_pFile,"P5\n%u\t%u\n%d\n",unsigned(m_ulWidth),unsigned(m_ulHeight),255);
  }

  if (m_ulY < m_ulHeight) {
    uint32_t x;
    for(x = 0;x < m_ulWidth;x++) {
      float v = line->Get(x) * m_dScale + m_dOffset;
      uint8_t b;

      if (v < 0) {
	b = 0;
      } else if (v > 255) {
	b = 255;
      } else {
	b = uint8_t(v);
      }

      fputc(b,m_pFile);
    }
    m_ulY++;
  }
}
///
