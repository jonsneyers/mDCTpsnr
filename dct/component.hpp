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
** $Id: component.hpp,v 1.7 2025/07/08 14:35:22 thor Exp $
**
*/

#ifndef DCT_COMPONENT_HPP
#define DCT_COMPONENT_HPP

/// Includes
#include <cstdint>
#include "global/thread.hpp"
#include <cassert>
///

/// Forwards
class Line;
///

/// class Component
// This class represents a single component of the image, line by line,
// in a overcomplete DCT basis.
class Component : private Thread {
  //
  // All lines we need for a DCT analysis, eight in total.
  class Line          *m_pInputLines[8];
  //
  // Output lines. Each frequency band gets one.
  class Line          *m_pOutputLines[8][8];
  //
  // Contiguous buffer for all 64 output lines (better cache locality)
  float              *m_pOutputBuffer;
  uint32_t                m_ulOutputStride;  // Stride between lines (aligned)
  //
  // Current Y position that is to be pushed in here.
  uint32_t                m_ulY;
  //
  // The window function.
  float               m_Window[8][8];
  //
  // The number of cores to use. 1 by default.
  int                  m_iCores;
  //
  // The original, unmodified line, offset to fit to the window.
  class Line          *m_pOriginal;
  //
  // The real DCT version.
  virtual void Run(int offset,int mod);
  //
  // Transform the buffered input lines into 64 bands using a lapped
  // window DCT.
  void DCT(void);
  //
public:
  Component(void);
  ~Component(void);
  //
  // Eagerly allocate all buffers once width is known (before any processing).
  // After this, all pointers are valid (never nullptr). Called before warmup loop.
  void AllocateBuffers(uint32_t width);
  //
  // Push a new line into the DCT conversion. Note that this is an overcomplete
  // representation, i.e. each line triggers a new DCT transformation.
  // PRECONDITION: AllocateBuffers() must be called first.
  void PushLine(Line &line);
  //
  // Get the current line indicating the given DCT band.
  // Returns nullptr only during warmup phase (first 7 image lines).
  // After warmup (line 8+), always returns valid Line* (never nullptr).
  // PRECONDITION: AllocateBuffers() must be called before any PushLine().
  class Line *GetDCTBand(int h,int v) const
  {
    assert(h >= 0 && h < 8 && v >= 0 && v < 8);

    return m_pOutputLines[h][v];
  }
  //
  // Return the offset-shifted original
  class Line *GetOriginal(void) const
  {
    return m_pOriginal;
  }
  //
  // Change the number of cores to use.
  // Must be set before the DCT is run.
  void SetCores(int cores)
  {
    m_iCores = cores;
  }
};
///

///
#endif
