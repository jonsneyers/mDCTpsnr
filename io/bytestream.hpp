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
** $Id: bytestream.hpp,v 1.2 2025/07/08 14:35:22 thor Exp $
**
*/

#ifndef IO_BYTESTREAM_HPP
#define IO_BYTESTREAM_HPP

/// Includes
#include <cassert>
#include <cstdint>
///

/// abstract ByteStream class
//
// The ByteStream class is an abstract data class to read or write groups
// or blocks of bytes from an abstract stream, additionally providing some
// kind of buffering mechanism. It implements non-virtual member functions
// that operate on the buffer only, and some virtual member functions to
// read or write complete buffers.
class ByteStream {
protected:
  uint32_t       m_ulBufSize;   // Size of our (internal) IO buffer
  uint32_t       m_ulBufBytes;  // # of valid(r)/available(w) bytes (0..bufsize-1)
  uint8_t      *m_pucBuffer;   // an IO buffer if we have it
  uint8_t      *m_pucBufPtr;   // a pointer to the first valid buffer byte
  uint32_t       m_ulCounter;   // counts output bytes, if possible
  //
  // Note: The counter and the buffers must be maintained by instances of
  // this abstract class.
  //
  // constructors: This just fills in the buffer size and resets the pointers
  //
  ByteStream(uint32_t bufsize = 2048)
    : m_ulBufSize(bufsize), m_ulBufBytes(0), m_pucBuffer(nullptr), m_ulCounter(0)
  { };
  //
  //
  //
  // fill up the buffer and flush it.
  // these two have to be replaced by the corresponding
  // member functions of the inherited classses
  virtual int32_t Fill(void) = 0;
  //
public: 
  // Seek modes for extended streams (though not this stream)
  enum SeekMode {
    Offset_Beginning = -1, // seek relative to the beginning
    Offset_Current   =  0, // seek relative to the current position
    Offset_End       =  1  // seek relative to the end of file.
  };
  //
  // The EOF indicator for get.
  static const int32_t Eof = -1;
  //
  // the destructor is virtual and has to be
  // overloaded.
  virtual ~ByteStream(void)
  { }
  //
  // Some rather standard IO functions, you know what they do.
  // These are not for overloading and non-virtual. All they
  // need is the buffer structure above. 
  int32_t Read(uint8_t *buffer,uint32_t size);          // read from buffer
  int32_t Write(const uint8_t *buffer,uint32_t size);   // write to buffer
  //
  // Reset the byte counter. This *MUST* be matched by a flush or a refill
  // or otherwise the result is undesirable.
  void ResetCounter(void)
  {
    m_ulCounter = 0;
  }
  // Flush the IO buffer. This must be defined by the instances of
  // this class.
  virtual void Flush(void) = 0;
  //
  // The following two methods are single byte IO functions,
  // inlined for maximal performance.
  int32_t Get(void)                          // read a single byte (inlined)
  {
    if (m_ulBufBytes == 0) {
      if (Fill() == 0)                    // Found EOF
	return Eof;
    }
    assert(int32_t(m_ulBufBytes) > 0);
    m_ulBufBytes--;
    return *m_pucBufPtr++;
  }
  //
  // Just the same for writing data.
  void Put(uint8_t byte)           // write a single byte (inlined)
  {
    if (m_ulBufBytes == 0) {
      Flush();                   // note that this will also allocate a buffer
    }
    assert(int32_t(m_ulBufBytes) > 0);
    m_ulBufBytes--;
    *m_pucBufPtr++ = byte;
  }
  //
  uint8_t LastByte(void)
  {
    if (m_pucBufPtr>m_pucBuffer) {
      return m_pucBufPtr[-1];
    }
    assert(0);
    return 0; // shut up
  }
  //
  // Return the last byte written/read and un-put/get it.
  uint8_t LastUnDo(void)
  {
    if (m_pucBufPtr>m_pucBuffer) {
      m_ulBufBytes++;
      m_pucBufPtr--;
      return *m_pucBufPtr;
    }
    assert(0);
    return 0; // shut-up g++
  }
  //
  // Return the byte counter = #of bytes read or written
  uint32_t FilePosition(void) 
  {
    if (m_pucBuffer) {
      return m_ulCounter + m_pucBufPtr - m_pucBuffer;
    } else {
      return m_ulCounter;
    }
  }
  //
};
///

///
#endif
