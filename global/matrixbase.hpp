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
** $Id: matrixbase.hpp,v 1.2 2025/07/08 14:35:22 thor Exp $
**
*/

#ifndef GLOBAL_MATRIXBASE_HPP
#define GLOBAL_MATRIXBASE_HPP

/// Includes
#include <cstdint>
#include "global/exceptions.hpp"
#include <cstdlib>
#include <cassert>
///

/// Matrix base class
class MatrixBase {
  //
protected:  
  // Memory management for the entry data. This holds the memory
  // to be kept here.
  struct MemKeeper {
    void*   m_pMem;       // where's the memory held?
    uint32_t  m_ulRefCount; // reference counter
    //
    MemKeeper(size_t size)
      : m_pMem(malloc(size)), m_ulRefCount(1)
    { 
      if (m_pMem == nullptr)
	Throw(NoMem,"MemKeeper::MemKeeper","out of memory");
    }
    //
    ~MemKeeper(void)
    {
      free(m_pMem);
    }
  }   *m_pMemory;
  //
  // This defines the number of entries I've to add to the index of
  // the (one dimensional!) data array to move from one row to the 
  // row below.
  // Note that this is in the coefficient type, so you'll need to multiply
  // by sizeof(type) to get the BytesPerRow.
  // Example:  A 4x4 (long) matrix usually has an EntriesPerRow==4.
  uint32_t m_ulEntriesPerRow;
  // The dimensions of the array. Humm, do we need to keep these?
  // At least for consistency checking, we should.
  uint32_t m_ulWidth,m_ulHeight;
  //
  // Build a new matrix allocating the indicated amount of storage.
  MatrixBase(size_t size)
    : m_pMemory(new MemKeeper(size))
  { 
    m_ulWidth         = 0;
    m_ulHeight        = 0;
    m_ulEntriesPerRow = 0;
  }
  //
  // Construct an empty matrix.
  MatrixBase(void)
    : m_pMemory(nullptr)
  { 
    m_ulWidth         = 0;
    m_ulHeight        = 0;
    m_ulEntriesPerRow = 0;
  }
  //
  // Copy construct a matrix from another source, reusing its
  // memory.
  MatrixBase(const MatrixBase &o)
  {
    if (o.m_pMemory) {
      __sync_add_and_fetch(&m_pMemory->m_ulRefCount,1); // used once more
      //m_pMemory->m_ulRefCount++; // used once more
      m_pMemory         = o.m_pMemory;
    } else {
      m_pMemory         = nullptr;
    }
    m_ulEntriesPerRow   = o.m_ulEntriesPerRow;
    m_ulWidth           = o.m_ulWidth;
    m_ulHeight          = o.m_ulHeight;
  }
  //
  // Assign the memory from another matrix to this matrix.
  MatrixBase &operator=(const MatrixBase &o)
  {
    // First increment the ref-count of the original
    // if that has memory. Otherwise, copy an empty
    // matrix.
    if (o.m_pMemory) {
      // This also protects for self-assignments.
      __sync_add_and_fetch(&o.m_pMemory->m_ulRefCount,1);
      //o.m_pMemory->m_ulRefCount++;
    }
    //
    // Now release our memory.
    if (m_pMemory && __sync_sub_and_fetch(&m_pMemory->m_ulRefCount,1) == 0) {
      delete m_pMemory;
    }
    m_pMemory         = o.m_pMemory;
    m_ulEntriesPerRow = o.m_ulEntriesPerRow;
    m_ulWidth         = o.m_ulWidth;
    m_ulHeight        = o.m_ulHeight;
    return *this;
  }
  //
  // Equip a matrix with its own memory
  void Allocate(size_t s)
  {
    // Release the current memory.
    if (m_pMemory && __sync_sub_and_fetch(&m_pMemory->m_ulRefCount,1) == 0) {
      delete m_pMemory; m_pMemory = nullptr;
    }
    m_pMemory         = new MemKeeper(s);
    m_ulWidth         = 0;
    m_ulHeight        = 0;
    m_ulEntriesPerRow = 0;
  }
  //
  // Release the contents of a matrix, making it empty again.
  void Release(void)
  {
    if (m_pMemory && __sync_sub_and_fetch(&m_pMemory->m_ulRefCount,1) == 0) {
      delete m_pMemory; m_pMemory = nullptr;
    }
    m_ulWidth         = 0;
    m_ulHeight        = 0;
    m_ulEntriesPerRow = 0;
  }
  //
  // Give up a matrix, release the memory reference.
  ~MatrixBase(void)
  {
    if (m_pMemory && __sync_sub_and_fetch(&m_pMemory->m_ulRefCount,1) == 0) {
      delete m_pMemory;
    }
  }
  //
public:
  // Return the coordinates given the index.
  // This should be avoided in frequent calling since it
  // requires multiplication and division.
  void CoordinatesOf(uint32_t idx,uint32_t &x,uint32_t &y) const
  {
    y = idx/m_ulEntriesPerRow;
    x = idx-y*m_ulEntriesPerRow;    
    assert(x < m_ulWidth && y < m_ulHeight);
  }
  //
  // The following method returns FALSE if the matrix is not yet
  // established
  bool IsEmpty(void) const
  {
    return (m_pMemory == nullptr);
  }
  //
  //
  // Simple dimension querry functions
  //
  uint32_t WidthOf(void) const
  {
    return m_ulWidth;
  }
  //
  uint32_t HeightOf(void) const
  {
    return m_ulHeight;
  }  
};
///

///
#endif
