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
#include <cstdint>
#include <cassert>
#include <cstring>
#include <stdlib.h>  // for aligned_alloc, free
///

/// class Line
// This class describes one horizontal line of the image data
class Line {
  //
  // The pointer to the allocated memory.
  float *m_pData;
  //
  // The nominal number of pixels, not inluding the boundary of the line.
  uint32_t   m_ulSize;
  //
  // Whether this Line owns its data buffer
  bool m_bOwnsData;
  //
  //
public:
  // Create a new line of the given nominal size.
  // Use 64-byte alignment for better cache/SIMD performance
  Line(uint32_t length)
    : m_pData(static_cast<float*>(aligned_alloc(64, ((length * sizeof(float) + 63) / 64) * 64)))
    , m_ulSize(length)
    , m_bOwnsData(true)
  { 
    if (!m_pData) m_pData = new float[length]; // fallback
  }
  //
  // Create a line that points to an external buffer (doesn't own it)
  Line(float *externalBuffer, uint32_t length)
    : m_pData(externalBuffer)
    , m_ulSize(length)
    , m_bOwnsData(false)
  { }
  //
  ~Line()
  {
    if (m_bOwnsData) {
      free(m_pData); // Use free() for aligned_alloc
    }
  }
  //
  // Copy from a second line.
  Line(const Line &org)
    : m_pData(static_cast<float*>(aligned_alloc(64, ((org.m_ulSize * sizeof(float) + 63) / 64) * 64)))
    , m_ulSize(org.m_ulSize)
    , m_bOwnsData(true)
  {
    if (!m_pData) m_pData = new float[org.m_ulSize];
    memcpy(m_pData,org.m_pData,m_ulSize * sizeof(float));
  }
  //
  // Assign from a second line.
  const Line &operator=(const Line &org)
  {
    float *tmp = static_cast<float*>(aligned_alloc(64, ((org.m_ulSize * sizeof(float) + 63) / 64) * 64));
    if (!tmp) tmp = new float[org.m_ulSize];
    memcpy(tmp,org.m_pData,org.m_ulSize * sizeof(float));
    if (m_bOwnsData) {
      free(m_pData); // Use free() for aligned_alloc
    }
    m_pData  = tmp;
    m_ulSize = org.m_ulSize;
    m_bOwnsData = true; // New allocation, we own it

    return *this;
  }
  //
  // Swap with the origin line.
  void Swap(Line &org)
  {
    float *tmp  = m_pData;
    uint32_t   siz  = m_ulSize;
    bool    own  = m_bOwnsData;

    m_ulSize     = org.m_ulSize;
    m_pData      = org.m_pData;
    m_bOwnsData  = org.m_bOwnsData;
    org.m_ulSize = siz;
    org.m_pData  = tmp;
    org.m_bOwnsData = own;
  }
  //
  // Copy the data from a second line in here, using an offset
  // in the source.
  void CopyFrom(const Line &org,size_t offset = 0)
  {
    assert(m_ulSize <= offset + org.m_ulSize);

    memcpy(m_pData,org.m_pData + offset,m_ulSize * sizeof(float));
  }
  //
  // Reset the line contents
  void Zero(void)
  {
    memset(m_pData,0,m_ulSize * sizeof(float));
  }
  //
  // Swap the data of this line with a different line
  //
  // Get the first nominal pixel of the line.
  float *Origin(void)
  {
    return m_pData;
  }
  //
  const float *Origin(void) const
  {
    return m_pData;
  }
  //
  // Return the pixel at the given offset.
  float &At(int offset)
  {
    assert(offset >= 0 && uint32_t(offset) < m_ulSize);

    return m_pData[offset];
  }
  //
  float Get(int offset) const
  { 
    assert(offset >= 0 && uint32_t(offset) < m_ulSize);

    return m_pData[offset];
  }
  //
  // Define the pixel at the given position.
  void Put(int offset,float v)
  {
    assert(offset >= 0 && uint32_t(offset) < m_ulSize);

    m_pData[offset] = v;
  }
  //
  // Return the size of the line in pixels, not including the extend.
  uint32_t LengthOf(void) const
  {
    return m_ulSize;
  }
};
///

///
#endif
