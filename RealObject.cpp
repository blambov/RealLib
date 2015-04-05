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

#include "defs.h"

#include <stdlib.h>
#include <limits.h>
#include <assert.h>
#include <string.h>
#include <math.h>

// STL
#include <vector>
#include <algorithm>

#include "RealObject.h"
#include "RealFuncs.h"

#ifndef NULL
#define NULL (0)
#endif

namespace RealLib {

ObjectList *g_pEstimatesList = NULL;


RealObject::RealObject(u32 rc)
: m_RefCount(rc), m_pPtrInObjList(0), m_pEstimate(), m_Depth(0), m_EstimateRefs(0)
{
}

RealObject::~RealObject()
{
    DestroyEstimate(); 
}

void RealObject::ReleaseSiblings()
{
}

RealObject* RealObject::AddRef()
{
    assert(this); 
    ++m_RefCount;

    // if there is an Encapsulation, the new reference
    // will need one more access to it
    if (m_EstimateRefs) 
        ++m_EstimateRefs;

    //    cout << "AddRef in " << this << " refs " << m_RefCount << " estrefs " << m_EstimateRefs << endl;
    return this; 
}

void RealObject::Release(int ReleaseCachedRefCount)
{
    assert(this); 
    assert(m_RefCount > 0);

    // note: we don't decrease m_EstimateRefs
    // this is because destruction happens
    // after the reference has received its evaluation
    if (m_EstimateRefs)
        m_EstimateRefs -= ReleaseCachedRefCount;
    if (m_EstimateRefs <= 0)
        DestroyEstimate();

    //    cout << "Release in " << this << " rem refs " << m_RefCount-1 << " estrefs " << m_EstimateRefs << endl;

    if (--m_RefCount == 0) {
        ReleaseSiblings();
        delete this;
    }
}

// GetEstimate(): returns an estimation of the term
// if it is referenced more than once, a copy of the
// evaluation is saved for later use

Encapsulation RealObject::GetEstimate()
{
    assert(this);
    assert(m_RefCount > 0);

    // is there a condition to destroy the Encapsulation? this is done by FinalizeRealLib
    //if (m_Precision != 0 && m_Precision < prec) DestroyEstimate();

    //    cout << "GetEstimate in " << this << " refs " << m_RefCount << " rem estrefs " << m_EstimateRefs-1 << endl;

    if (m_pEstimate) {
        // if this is the last time we'll need the Encapsulation
        if (--m_EstimateRefs == 0) {
            Encapsulation est = *m_pEstimate;
            DestroyEstimate();
            return est;
        } else return *m_pEstimate;
    } else if (m_RefCount > 1) {
        // create a new Encapsulation
        Encapsulation val(Evaluate());
        m_pEstimate = CreateEstimate(val);
        m_EstimateRefs = m_RefCount - 1;

        // return it
        return val;
    } else {
        // this object is referenced only once
        // don't save, but clear m_EstimateRefs
        // which may have been set by HasEstimateReference
        m_EstimateRefs = 0;
        return Evaluate();
    }

}

void RealObject::AddToEstimatesList()
{
    // save object in g_pEstimatesList
    m_pPtrInObjList = new ObjectList;
    assert(m_pPtrInObjList);
    m_pPtrInObjList->next = g_pEstimatesList;
    m_pPtrInObjList->prev = NULL;
    m_pPtrInObjList->obj = this;
    if (g_pEstimatesList)
        g_pEstimatesList->prev = m_pPtrInObjList;
    g_pEstimatesList = m_pPtrInObjList;
}

void RealObject::RemoveFromEstimatesList()
{
    // delete it from g_pEstimatesList
    assert(m_pPtrInObjList);
    if (m_pPtrInObjList->prev) {
        m_pPtrInObjList->prev->next = m_pPtrInObjList->next;
    } else {
        g_pEstimatesList = m_pPtrInObjList->next;
    }
    if (m_pPtrInObjList->next) {
        m_pPtrInObjList->next->prev = m_pPtrInObjList->prev;
    }
    delete m_pPtrInObjList;

    m_pPtrInObjList = NULL;
}

EncapsulationPointer RealObject::CreateEstimate(const Encapsulation &val)
{
    AddToEstimatesList();
    return EncapsulationPointer(val);
}

// DestroyEstimate(): if the value is cached, destroy it
// so reevaluation can begin

void RealObject::DestroyEstimate()
{
    if (m_pEstimate) {
        //        cout << "destroying estimate " << *m_pEstimate << endl;
        m_pEstimate.Release();

        RemoveFromEstimatesList();
        m_EstimateRefs = 0;
    }
}

RealObject* RealObject::GetSibling(int index)
{
    index;
    return NULL;
}

void RealObject::NonRecursiveRelease()
{
    // release the object, but do not release siblings if it
    // is deleted. They have already been visited 
    // (FinishRelease is called on GetDepthList exit,
    // from the bottom up)
    assert(this);
    assert(m_RefCount >= 1);

    if (--m_RefCount == 0)
        delete this;
}

// real from double

RealFromDouble::RealFromDouble(const double src)
: m_Value(src)
{
}

RealFromDouble::~RealFromDouble()
{
}

Encapsulation RealFromDouble::Evaluate()
{
    return Encapsulation(m_Value);
}

// real from string

RealFromString::RealFromString(const char *src)
{
    // the source needs to be copied
    // otherwise the user might change it later
    assert(src);
    m_pString = new char[strlen(src) + 1];
    assert(m_pString);
#ifdef _MSC_VER
#pragma warning (push)
#pragma warning (disable:4996)
    strcpy(m_pString, src);
#pragma warning (pop)
#else
    strcpy(m_pString, src);
#endif
}

RealFromString::~RealFromString()
{
    delete m_pString;
}

Encapsulation RealFromString::Evaluate()
{
    return Encapsulation(m_pString);
}

// real from oracle function

RealFromOracle::RealFromOracle(OracleFunction oracle)
: m_pOracle(oracle)
{
}

RealFromOracle::~RealFromOracle()
{
}

Encapsulation RealFromOracle::Evaluate()
{
    const char *val = m_pOracle(g_WorkingPrecision);
    const char *eptr = strchr(val, 'e');
    if (!eptr) eptr = strchr(val, 'E');
    if (!eptr) eptr = strchr(val, 0);   // strlen of a kind
    i32 len;
    if (strchr(val, '.')) len = i32(eptr - val - 1);
    else len = i32(eptr - val);
    double exp = -len * LOG_2_10;
    double man = (exp - floor(exp)) / LOG_2_10;
    Encapsulation e(val);
    Encapsulation err(man);
    err = err << i32(floor(exp));
    return e.AddError(err);
}

// real nullary (constants)

RealNullary::RealNullary(FuncNullary pFunc, UserInt user)
: m_pFunc(pFunc), m_iUserData(user)
{
    assert(pFunc);
}

RealNullary::~RealNullary()
{
}

Encapsulation RealNullary::Evaluate()
{
    return m_pFunc(g_WorkingPrecision, m_iUserData);
}

// real unary

RealUnary::RealUnary(FuncUnary pFunc, RealObject *pArg, UserInt user)
: m_pFunc(pFunc), m_pArg(pArg), m_iUserData(user)
{
    assert(pFunc);
    assert(pArg);

    // reference the sibling
    pArg->AddRef();

    // mark the depth
    m_Depth = pArg->m_Depth + 1;
}

RealUnary::~RealUnary()
{
}

void RealUnary::ReleaseSiblings()
{
    m_pArg->Release();
}

// evaluate: evaluate argument, apply function to it
Encapsulation RealUnary::Evaluate()
{
    return m_pFunc(m_pArg->GetEstimate(), m_iUserData);
}

RealObject* RealUnary::GetSibling(int index)
{
    return index == 0 ? m_pArg : NULL;
}


// binary

RealBinary::RealBinary(FuncBinary pFunc, RealObject *pLeft, RealObject *pRight, UserInt user)
: m_pFunc(pFunc), m_pLeft(pLeft), m_pRight(pRight), m_iUserData(user)
{
    assert(pFunc);
    assert(pLeft);
    assert(pRight);

    // reference the siblings
    pLeft->AddRef();
    pRight->AddRef();

    // set depth (greater than both siblings)
    m_Depth = max(pLeft->m_Depth, pRight->m_Depth) + 1;
}

RealBinary::~RealBinary()
{
}

void RealBinary::ReleaseSiblings()
{
    m_pLeft->Release();
    m_pRight->Release();
}

Encapsulation RealBinary::Evaluate()
{
    return m_pFunc(m_pLeft->GetEstimate(), m_pRight->GetEstimate(), m_iUserData);
}

RealObject* RealBinary::GetSibling(int index)
{
    if (index == 0) {
        if (m_pLeft->m_Depth <= m_pRight->m_Depth) return m_pRight;
        else return m_pLeft;
    } else if (index == 1) {
        if (m_pLeft->m_Depth <= m_pRight->m_Depth) return m_pLeft;
        else return m_pRight;
    } else return NULL;

    //return index == 0 ? m_pLeft : (index == 1 ? m_pRight : NULL);
}

// RealArray implementation

RealArray::RealArray(FuncArray pFunc, RealObject **pArray, unsigned int count, UserInt userdata)
: m_pFunc(pFunc), m_pArray(pArray), m_uCount(count), m_iUserData(userdata)
{
    assert(pFunc);
    assert(pArray);
    // references should have already been added

    u32 depth = 0;
    for (u32 i=0;i<count;++i)
        depth = max(i32(depth), pArray[i]->m_Depth);
    m_Depth = depth + 1;
}

RealArray::~RealArray()
{
    DestroyEstimate();
    //   assert(!m_pEstimateArray);
    delete m_pArray;
}

void RealArray::SetRequestIndex(i32 index)
{
    assert(index >= 0 && index < i32(m_uCount)); 
    m_uRequestIndex = index; 
    if (m_pEstimate) m_pEstimate = m_pEstimateArray + index; 
}

void RealArray::ReleaseSiblings()
{
    for (u32 i=0;i<m_uCount;++i)
        m_pArray[i]->Release();
}

EncapsulationPointer RealArray::CreateEstimate(const Encapsulation &val)
{
    val;
    // the actual job is done by Evaluate
    // only return the pointer
    return EncapsulationPointer::FromPointer
            (m_pEstimateArray + m_uRequestIndex);
}

void RealArray::DestroyEstimate()
{
    if (m_pEstimateArray) {
        m_pEstimateArray.ReleaseArray(m_uCount);
        m_pEstimate = EncapsulationPointer();
        m_EstimateRefs = 0;
        RemoveFromEstimatesList();
    }
}

Encapsulation RealArray::Evaluate()
{
    m_pEstimateArray = EncapsulationPointer(m_uCount);
    AddToEstimatesList();

    for (u32 i=0;i<m_uCount;++i)
        m_pEstimateArray[i] = m_pArray[i]->GetEstimate();
    ArrayInterface<Encapsulation> arr(m_pEstimateArray, m_uCount);
    m_pFunc(arr, m_iUserData);

    return m_pEstimateArray[m_uRequestIndex];
}

RealObject* RealArray::GetSibling(int index)
{
    return u32(index) < m_uCount ? m_pArray[index] : NULL;
}

RealArrayElement::RealArrayElement(RealArray *pArray, u32 myindex)
: m_pArray(pArray), m_uMyIndex(myindex)
{
    pArray->AddRef();
    m_Depth = (pArray->m_Depth + 1);
}

RealArrayElement::~RealArrayElement()
{
}

void RealArrayElement::ReleaseSiblings()
{
    m_pArray->Release();
}

EncapsulationPointer RealArrayElement::CreateEstimate(const Encapsulation &val)
{
    // the Encapsulation is already saved in the array
    return m_pArray->CreateEstimate(val);
}

void RealArrayElement::DestroyEstimate()
{
    // do nothing, the array will clean it
    m_pEstimate = EncapsulationPointer();
}

Encapsulation RealArrayElement::Evaluate()
{
    m_pArray->SetRequestIndex(m_uMyIndex);
    return m_pArray->GetEstimate();
}

RealObject* RealArrayElement::GetSibling(int index)
{
    return index == 0 ? m_pArray : NULL;
}

}  // namespace
