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
** $Id: line.hpp,v 1.6 2025/07/08 14:35:22 thor Exp $
**
*/

#ifndef DCT_LINE_HPP
#define DCT_LINE_HPP

/// Includes
#include "global/types.hpp"
#include "std/assert.hpp"
#include "std/string.hpp"
///

/// class Line
// This class describes one horizontal line of the image data
class Line {
  //
  // The pointer to the allocated memory.
  DOUBLE *m_pData;
  //
  // The nominal number of pixels, not inluding the boundary of the line.
  ULONG   m_ulSize;
  //
  //
public:
  // Create a new line of the given nominal size.
  Line(ULONG length)
    : m_pData(new DOUBLE[length]), m_ulSize(length)
  { }
  //
  ~Line()
  {
    delete[] m_pData;
  }
  //
  // Copy from a second line.
  Line(const Line &org)
    : m_pData(new DOUBLE[org.m_ulSize]), m_ulSize(org.m_ulSize)
  {
    memcpy(m_pData,org.m_pData,m_ulSize * sizeof(DOUBLE));
  }
  //
  // Assign from a second line.
  const Line &operator=(const Line &org)
  {
    DOUBLE *tmp = new DOUBLE[org.m_ulSize];
    memcpy(tmp,org.m_pData,org.m_ulSize * sizeof(DOUBLE));
    delete[] m_pData;
    m_pData  = tmp;
    m_ulSize = org.m_ulSize;

    return *this;
  }
  //
  // Swap with the origin line.
  void Swap(Line &org)
  {
    DOUBLE *tmp  = m_pData;
    ULONG   siz  = m_ulSize;

    m_ulSize     = org.m_ulSize;
    m_pData      = org.m_pData;
    org.m_ulSize = siz;
    org.m_pData  = tmp;
  }
  //
  // Copy the data from a second line in here, using an offset
  // in the source.
  void CopyFrom(const Line &org,size_t offset = 0)
  {
    assert(m_ulSize <= offset + org.m_ulSize);

    memcpy(m_pData,org.m_pData + offset,m_ulSize * sizeof(DOUBLE));
  }
  //
  // Reset the line contents
  void Zero(void)
  {
    memset(m_pData,0,m_ulSize * sizeof(DOUBLE));
  }
  //
  // Swap the data of this line with a different line
  //
  // Get the first nominal pixel of the line.
  DOUBLE *Origin(void)
  {
    return m_pData;
  }
  //
  const DOUBLE *Origin(void) const
  {
    return m_pData;
  }
  //
  // Return the pixel at the given offset.
  DOUBLE &At(int offset)
  {
    assert(offset >= 0 && ULONG(offset) < m_ulSize);

    return m_pData[offset];
  }
  //
  DOUBLE Get(int offset) const
  { 
    assert(offset >= 0 && ULONG(offset) < m_ulSize);

    return m_pData[offset];
  }
  //
  // Define the pixel at the given position.
  void Put(int offset,DOUBLE v)
  {
    assert(offset >= 0 && ULONG(offset) < m_ulSize);

    m_pData[offset] = v;
  }
  //
  // Return the size of the line in pixels, not including the extend.
  ULONG LengthOf(void) const
  {
    return m_ulSize;
  }
};
///

///
#endif
