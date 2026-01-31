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
** $Id: filestream.hpp,v 1.2 2025/07/08 14:35:22 thor Exp $
**
*/

#ifndef IO_FILESTREAM_HPP
#define IO_FILESTREAM_HPP

/// Includes
#include "io/bytestream.hpp"
#include <cstdio>
///

/// FileStream
// This class provides binary file reading from
// the Os environment, or writing to it.
class FileStream : public ByteStream {
  //
  // The stdio file pointer we use for reading the
  // file.
  FILE *m_pFile;
  //
  // Set if the file is open for writing.
  bool  m_bOpenForWrite;
  //
  // Methods required for the bytestream implementation.
  // these two have to be replaced by the corresponding
  // member functions of the inherited classses
  virtual int32_t Fill(void);
  //  
  // Flush the IO buffer.
  virtual void Flush(void);
  //
public:
  FileStream(uint32_t bufsize = 2048)
    : ByteStream(bufsize), m_pFile(nullptr)
  { }
  //
  virtual ~FileStream(void);
  //
  // Open a stream for input. Throws on error.
  void OpenForRead(const char *path,uint32_t bufsize = 2048);
  //
  // Open a stream for writing. Throws on error.
  void OpenForWrite(const char *path,uint32_t bufsize = 2048);
  //
  // Close a file (also happens on destruction)
  void Close(void);
};
///

///
#endif
