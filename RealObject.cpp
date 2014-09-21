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
: m_RefCount(rc), m_pPtrInObjList(0), m_pEstimate(NULL), m_Depth(0), m_EstimateRefs(0)
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
    return this; 
}

void RealObject::Release() 
{ 
    assert(this); 
   assert(m_RefCount > 0);

    // note: we don't decrease m_EstimateRefs
    // this is because destruction happens
    // after the reference has received its evaluation

    if (--m_RefCount == 0) {
        ReleaseSiblings();
        delete this; 
    }
}

// Encapsulation(prec): returns an estimation of the term
// if it is referenced more than once, a copy of the
// evaluation is saved for later use

Encapsulation RealObject::GetEstimate()
{
    assert(this);
    assert(m_RefCount > 0);

    // is there a condition to destroy the Encapsulation? this is done by FinalizeRealLib
   //if (m_Precision != 0 && m_Precision < prec) DestroyEstimate();

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

Encapsulation *RealObject::CreateEstimate(const Encapsulation &val)
{
   AddToEstimatesList();
   return new Encapsulation(val);
}

// DestroyEstimate(): if the value is cached, destroy it
// so reevaluation can begin

void RealObject::DestroyEstimate()
{ 
    if (m_pEstimate) {
        delete m_pEstimate; 
        m_pEstimate = NULL;

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

/* obsolete
bool RealObject::HasEstimateReference(i32 Depth)
{
    // return true if there is an Encapsulation, or if the object
    // is already visited.
    if (m_EstimateRefs == 0 && m_Depth >= Depth) {
        m_EstimateRefs = 1;
        return false;
    } else return true;
}

bool RealObject::StartRelease(i32 Depth)
{
   if (m_Depth < Depth) {
      // release siblings so that combined with FinishRelease it does a full release
      if (--m_RefCount == 0) ReleaseSiblings();
      return true;
   }
   return --m_RefCount > 0;  
}

void RealObject::FinishRelease()
{
    // release the object, but do not release siblings if it
    // is deleted. They have already been visited 
    // (FinishRelease is called on GetDepthList exit,
    // from the bottom up)
    assert(this);
    assert(m_RefCount >= 0);

    if (m_RefCount == 0)
        delete this;
}

void RealObject::GetDepthList(i32 Depth, TesterFunc Test, OnExitFunc OnExit, EvalList **pList, bool DiveIntoDeeperFirst)
{
   pList; DiveIntoDeeperFirst;
    // no siblings. only take care to call Test and OnExit
    (this->*Test)(Depth);
   if (OnExit) (this->*OnExit)();
}
*/

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
    strcpy(m_pString, src);
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

/* obsolete
void RealUnary::GetDepthList(i32 Depth, TesterFunc Test, OnExitFunc OnExit, EvalList **pList, bool DiveIntoDeeperFirst)
{
    // do we need to continue?
    if (!(Test && (this->*Test)(Depth))) {
       // have we reached maximum depth?
       if (Depth <= 1) {
           if (pList) {
               EvalList *pNew = new EvalList;
               pNew->next = *pList;
               pNew->obj = AddRef();
               pNew->used = false;
               *pList = pNew;
           }
       } else {
           // continue with sibling
           m_pArg->GetDepthList(Depth-1, Test, OnExit, pList, DiveIntoDeeperFirst);
       }
   }

   // clean up
   if (OnExit) (this->*OnExit)();
}
*/
// binary on int
/*
RealBinaryOnInt::RealBinaryOnInt(FuncBinaryOnInt pFunc, RealObject *pArg, long iarg)
: m_pFunc(pFunc), m_pArg(pArg), m_iArg(iarg)
{
    assert(pFunc);
    assert(pArg);
    
    // reference the sibling
    pArg->AddRef();

    // mark the depth
    m_Depth = pArg->m_Depth + 1;
}

RealBinaryOnInt::~RealBinaryOnInt()
{
}

void RealBinaryOnInt::ReleaseSiblings()
{
    m_pArg->Release();
}

// evaluate: evaluate argument, apply function to it
Encapsulation RealBinaryOnInt::Evaluate()
{
    return m_pFunc(m_pArg->GetEstimate(), m_iArg);
}

RealObject* RealBinaryOnInt::GetSibling(int index)
{
    return index == 0 ? m_pArg : NULL;
}*/

/* obsolete
void RealBinaryOnInt::GetDepthList(i32 Depth, TesterFunc Test, OnExitFunc OnExit, EvalList **pList, bool DiveIntoDeeperFirst)
{
    // do we need to continue?
    if (!(Test && (this->*Test)(Depth))) {
       // have we reached maximum depth?
       if (Depth <= 1) {
           if (pList) {
               EvalList *pNew = new EvalList;
               pNew->next = *pList;
               pNew->obj = AddRef();
               pNew->used = false;
               *pList = pNew;
           }
       } else {
           // continue with sibling
           m_pArg->GetDepthList(Depth-1, Test, OnExit, pList, DiveIntoDeeperFirst);
       }
   }

   // clean up
   if (OnExit) (this->*OnExit)();
}
*/

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
    return index == 0 ? m_pLeft : index == 1 ? m_pRight : NULL;
}

/* obsolete
void RealBinary::GetDepthList(i32 Depth, TesterFunc Test, OnExitFunc OnExit, EvalList **pList, bool DiveIntoDeeperFirst)
{
    // do we need to continue?
    if (!(Test && (this->*Test)(Depth))) {
       // have we reached maximum depth?
       if (Depth <= 1) {
           if (pList) {
               EvalList *pNew = new EvalList;
               pNew->next = *pList;
               pNew->obj = AddRef();
               pNew->used = false;
               *pList = pNew;
           }
       } else {
           // continue with the sibling with greater or lower depth first
         // depending on DiveIntoDeeperFirst flag.
           // otherwise a previous run might have marked visited
           // something that could lead to a deeper dive on evaluation
           if ((m_pLeft->m_Depth >= m_pRight->m_Depth) == DiveIntoDeeperFirst) {
               m_pLeft->GetDepthList(Depth-1, Test, OnExit, pList, DiveIntoDeeperFirst);
               m_pRight->GetDepthList(Depth-1, Test, OnExit, pList, DiveIntoDeeperFirst);
           } else {
               m_pRight->GetDepthList(Depth-1, Test, OnExit, pList, DiveIntoDeeperFirst);
               m_pLeft->GetDepthList(Depth-1, Test, OnExit, pList, DiveIntoDeeperFirst);
           }
      }
   }

   // clean up
   if (OnExit) (this->*OnExit)();
}
*/

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

Encapsulation* RealArray::CreateEstimate(const Encapsulation &val)
{
   val;
   // the actual job is done by Evaluate
   // only return the pointer
   return m_pEstimateArray + m_uRequestIndex;
}

void RealArray::DestroyEstimate()
{
   if (m_pEstimateArray) {
      delete [] m_pEstimateArray;
      m_pEstimateArray = m_pEstimate = NULL;
        m_EstimateRefs = 0;
      RemoveFromEstimatesList();
   }
}

Encapsulation RealArray::Evaluate()
{
   m_pEstimateArray = new Encapsulation[m_uCount];
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

/* obsolete
bool RealObjectDeeper(RealObject *left, RealObject *right)
{
   return left->m_Depth >= right->m_Depth;
}

void RealArray::GetDepthList(i32 Depth, TesterFunc Test, OnExitFunc OnExit, EvalList **pList, bool DiveIntoDeeperFirst)
{
    // do we need to continue?
    if (!(Test && (this->*Test)(Depth))) {
       // have we reached maximum depth?
       if (Depth <= 1) {
           if (pList) {
               EvalList *pNew = new EvalList;
               pNew->next = *pList;
               pNew->obj = AddRef();
               pNew->used = false;
               *pList = pNew;
           }
       } else {
           // continue with siblings
         using namespace std;

//         vector<RealObject*> arr(vector<RealObject*>::iterator(m_pArray), 
//                                 vector<RealObject*>::iterator(m_pArray + m_uCount));//
                                 
         vector<RealObject*> arr(m_uCount);
         for (u32 i=0;i<m_uCount;++i)
            arr.push_back(m_pArray[i]);
         sort(arr.begin(), arr.end(), RealObjectDeeper);

         if (DiveIntoDeeperFirst)
            for (vector<RealObject*>::iterator it=arr.begin();it!=arr.end();++it)
                 (*it)->GetDepthList(Depth-1, Test, OnExit, pList, DiveIntoDeeperFirst);
         else
            for (vector<RealObject*>::reverse_iterator it=arr.rbegin();it!=arr.rend();++it)
                 (*it)->GetDepthList(Depth-1, Test, OnExit, pList, DiveIntoDeeperFirst);
       }
   }

   // clean up
   if (OnExit) (this->*OnExit)();
}
*/

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

Encapsulation* RealArrayElement::CreateEstimate(const Encapsulation &val)
{
   // the Encapsulation is already saved in the array
   return m_pArray->CreateEstimate(val);
}

void RealArrayElement::DestroyEstimate()
{
   // do nothing, the array will clean it
   m_pEstimate = NULL;
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

/* obsolete
void RealArrayElement::GetDepthList(i32 Depth, TesterFunc Test, OnExitFunc OnExit, EvalList **pList, bool DiveIntoDeeperFirst)
{
    // do we need to continue?
    if (!(Test && (this->*Test)(Depth))) {
       // have we reached maximum depth?
       if (Depth <= 1) {
           if (pList) {
               EvalList *pNew = new EvalList;
               pNew->next = *pList;
               pNew->obj = AddRef();
               pNew->used = false;
               *pList = pNew;
           }
       } else {
           // continue with sibling
           m_pArray->GetDepthList(Depth-1, Test, OnExit, pList, DiveIntoDeeperFirst);
       }
   }

   // clean up
   if (OnExit) (this->*OnExit)();
}
*/

/*
bool RealObject::HasEstimate()
{
    return m_pEstimate != 0;
}

void RealObject::SetEstimateReference()
{
    ++m_EstimateRefs;
}


*/

}  // namespace
