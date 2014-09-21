/*

  RealLib, a library for efficient exact real computation
  Copyright (C) 2006  Branimir Lambov <barnie@brics.dk>

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

*/ 

/*

  DataManager.h

  Wrapper around implementation details that handle the allocation of
  space for mantissas.

*/

#ifndef FILE_DATA_MANAGER_H
#define FILE_DATA_MANAGER_H

#include "defs.h"
#include <assert.h>

namespace RealLib {

// Alloc: object that handles multiple references to a memory block
typedef u32 Alloc;

// DataBuffer holds all the memory blocks. It can grow if necessary.

class DataBuffer {
private:
    u32 m_uCount;
    u32 m_uPrec;
    u32 *m_pData;

public:
    DataBuffer(u32 howmany, u32 precision)
        : m_uCount(howmany), m_uPrec(precision), m_pData(new u32 [howmany * precision])
    {}

    ~DataBuffer() { if (m_pData) delete [] m_pData; }

    bool isValid() { return !!m_pData; }

    bool Grow (u32 howmuch);
    u32 GetSize() { return m_uCount; }
    u32 GetPrec() { return m_uPrec; }

    // get a pointer to the words of the indexth block of memory
    u32* operator [] (u32 index) { assert(index < m_uCount); return m_pData + index * m_uPrec; }

};

// FreeStack: a stack that holds all available memory block indices.
// on initialization all indices from 0 to m_uSize-1 are enumerated.
// an index is popped on allocation and pushed on free.
class FreeStack {
private:
    u32 *m_pData;
    u32 m_uCount;
    u32 m_uSize;

public:
    FreeStack(u32 howmany);
    ~FreeStack();

    bool isValid() { return !!m_pData; }

    void push(u32 index) 
    {   assert(m_uCount < m_uSize); 
        m_pData[m_uCount++] = index; }
    u32 pop()
    {   assert(m_uCount > 0);
        return m_pData[--m_uCount]; }

    bool Grow(u32 howmuch);
    u32 GetCount() { return m_uCount; }
    u32 GetSize() { return m_uSize; }
};

// DataManager: an abstaction layer for data allocation and freeing
// uses DataBuffer and FreeStack
// can be implemented differently
class DataManager {
private:
    DataBuffer m_Buf;
    FreeStack m_Free;
    u32 m_uGrow;

    // allocate a new mantissa by popping an index from m_Free
    // grow if necessary
    u32 get();
    // free by pushing index to m_Free
    void free(u32 index);

public:
    // constructor arguments: precision, initial size, and growth step
    DataManager(u32 prec, u32 howmany = 100, u32 grow = 100);

    bool isValid() { return m_Buf.isValid() && m_Free.isValid(); }

    u32 GetFreeCount() { return m_Free.GetCount(); }
    u32 GetSize() { return m_Free.GetSize(); }
    u32 GetPrec() { return m_Buf.GetPrec(); }

    // allocation of new block
    Alloc newAlloc();
    // referencing a block
    Alloc referenceAlloc(Alloc alloc);
    // releasing it
    void releaseAlloc(Alloc alloc);

    // hide implementation details
    u32 *operator [] (Alloc index) { return m_Buf[index] + 1; }

    // allocs with more than one reference can't be changed
    bool AllocCanBeChanged (Alloc alloc) { return m_Buf[alloc][0] == 1; }
};

} // namespace

#endif // file
