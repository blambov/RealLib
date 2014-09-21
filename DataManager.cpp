/*

  RealLib, a library for efficient exact real computation
  Copyright (C) 2006 Branimir Lambov

  This library is licensed under the Apache License, Version 2.0 (the "License");
  you may not use this library except in compliance with the License.
  You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.
*/

#include <stdlib.h>

#include "DataManager.h"

namespace RealLib {

// grow to accomodate more data
bool DataBuffer::Grow(u32 howmuch)
{
    u32* pNew = new u32 [(m_uCount + howmuch) * m_uPrec];
    if (!pNew) return false;

    for (u32 u = 0; u < m_uCount * m_uPrec; ++u)
        pNew[u] = m_pData[u];
    delete [] m_pData;

    m_pData = pNew;
    m_uCount += howmuch;
    return true;
}

// initialize by adding all indices
FreeStack::FreeStack(u32 howmany)
: m_pData(new u32 [howmany]), m_uCount(howmany), m_uSize(howmany)
{
    if (m_pData)
        for (u32 u = 0; u<howmany; ++u) m_pData[u] = howmany-u-1;
}

FreeStack::~FreeStack()
{
    // only allow release if all data is freed
    assert(m_uCount == m_uSize);
    if (m_pData) delete [] m_pData;
}

// grow and add all new indices
bool FreeStack::Grow(u32 howmuch)
{
    assert(isValid());
    u32* pNew = new u32 [m_uSize + howmuch];
    if (!pNew) return false;

    u32 u;
    
    for (u = 0; u < howmuch; ++u)
        pNew[u] = m_uSize + howmuch - u - 1;
    
    // this wouldn't be usually needed 
    // (Grow should be invoked at m_uCount==0)
    // do it anyway
    for (; u < m_uCount; ++u)
        pNew[u] = m_pData[u - howmuch];

    delete [] m_pData;

    m_pData = pNew;
    m_uSize += howmuch;
    m_uCount += howmuch;
    return true;
}

// Data manager constructor
// allocates space for data and reference counts in m_Buf
// and for free stack in m_Free
DataManager::DataManager(u32 prec, u32 howmany, u32 grow)
: m_Buf(howmany, prec+1), m_Free(howmany), m_uGrow(grow)
{
}

u32 DataManager::get()
{
    if (GetFreeCount() == 0)
        if (!m_Buf.Grow(m_uGrow) || !m_Free.Grow(m_uGrow)) 
            return u32(-1);     // should be exception
        
    return m_Free.pop();
}

void DataManager::free(u32 index)
{
    m_Free.push(index);
}

Alloc DataManager::newAlloc()
{
    Alloc r = get();
    // the first word of the block will hold the reference count
    m_Buf[r][0] = 1;
    return r;
}

Alloc DataManager::referenceAlloc(Alloc alloc)
{
    ++m_Buf[alloc][0];
    return alloc;
}

void DataManager::releaseAlloc(Alloc alloc)
{
    assert(i32(m_Buf[alloc][0]) > 0);

    if (--m_Buf[alloc][0] == 0) free(alloc);
}

} // namespace
