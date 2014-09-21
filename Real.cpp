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
#include <string.h>
#include <assert.h>
#include <math.h>
#include <iostream>
#include <string>
#include <vector>
#include <valarray>
#include "Real.h"
#include "RealFuncs.h"
#include "RealObject.h"

namespace RealLib {

unsigned g_precMax = 0;
unsigned g_memReq = 0;

//const Real Pi(pi);
//const Real Ln2(ln2);
    
void InitializeRealLib(unsigned precStart, unsigned precMax, unsigned numEst)
{
    // set maximum precision
    g_precMax = precMax;
    g_memReq = precStart * numEst;

    // initialize the underlying arbitrary precision engine
    InitializeLongFloat(precStart, numEst);
}

unsigned FinalizeRealLib()
{
    unsigned res = g_WorkingPrecision;
    g_precMax = 0;

    // release all temporary estimates
    while (g_pEstimatesList) {
        g_pEstimatesList->obj->DestroyEstimate();
    }
   DestroyConstantEstimates();
    
    // finalize the arbitrary precision engine
    FinalizeLongFloat();

    return res;
}

unsigned ResetRealLib(unsigned precStart)
{
    unsigned precMax = g_precMax;
    unsigned res = FinalizeRealLib();
    InitializeRealLib(precStart, precMax, (g_memReq + precStart - 1) / precStart);
    return res;
}

Real::Real(const Real &src)
: m_pObject(src.m_pObject)
{
    assert(m_pObject);
    m_pObject->AddRef();
}

Real::Real(RealObject *src)
: m_pObject(src)
{
    // copy the object and increase its reference count
    assert(m_pObject);
    m_pObject->AddRef();
}

Real::Real(const char *src)
{
    // create a new RealFromString object and link to it
    m_pObject = new RealFromString(src);
    assert(m_pObject);
    m_pObject->AddRef();
}

Real::Real(const double src)
{
    // create a new RealFromDouble object and link to it
    m_pObject = new RealFromDouble(src);
    assert(m_pObject);
    m_pObject->AddRef();
}

Real::Real(OracleFunction oracle)
{
    m_pObject = new RealFromOracle(oracle);
    assert(m_pObject);
    m_pObject->AddRef();
}

Real::Real(Real::FuncNullary constant, UserInt user)
{
    // create a new RealNullary object and link to it
    m_pObject = new RealNullary(constant, user);
    assert(m_pObject);
    m_pObject->AddRef();
}
/*
Real::Real(Real::FuncUnary unfunc, const Real &arg)
{
    // create a new RealUnary object and link to it
   m_pObject = new RealUnary(unfunc, arg.m_pObject, 0);
    assert(m_pObject);
    m_pObject->AddRef();
}

Real::Real(Real::FuncBinary binfunc, const Real &lhs, const Real &rhs)
{
   // create a new RealBinary object and link to it
   m_pObject = new RealBinary(binfunc, lhs.m_pObject, rhs.m_pObject, 0);
    assert(m_pObject);
    m_pObject->AddRef();
}*/

Real::Real(Real::FuncUnary unfunc, const Real &arg, UserInt user)
{
    // create a new RealUnary object and link to it
   m_pObject = new RealUnary(unfunc, arg.m_pObject, user);
    assert(m_pObject);
    m_pObject->AddRef();
}

Real::Real(Real::FuncBinary binfunc, const Real &lhs, const Real &rhs, UserInt user)
{
   // create a new RealBinary object and link to it
   m_pObject = new RealBinary(binfunc, lhs.m_pObject, rhs.m_pObject, user);
    assert(m_pObject);
    m_pObject->AddRef();
}

/*
Real::Real(Real::FuncBinaryOnInt binfunc, const Real &lhs, long rhs)
{
   // create a new RealBinaryOnInt object and link to it
    m_pObject = new RealBinaryOnInt(binfunc, lhs.m_pObject, rhs);
    assert(m_pObject);
    m_pObject->AddRef();
}*/

template <class ARRAY>
void RealApplyArrayFunction(RealFuncArray arfunc, ARRAY &arr, UserInt user)
{
    long count = (long)arr.size();
   RealObject **parr = new RealObject* [count];
   for (long i=0;i<count;++i) {
      parr[i] = arr[i].m_pObject;
      parr[i]->AddRef();
   }
   RealArray *oarr = new RealArray(arfunc, parr, count, user);
   assert(oarr);

   for (long i=0;i<count;++i)
      arr[i] = Real(new RealArrayElement(oarr, i));
}

template
void RealApplyArrayFunction(RealFuncArray arfunc, ArrayInterface<Real> &arr, UserInt userdata);

template
void RealApplyArrayFunction(RealFuncArray arfunc, std::valarray<Real> &arr, UserInt userdata);

template
void RealApplyArrayFunction(RealFuncArray arfunc, std::vector<Real> &arr, UserInt userdata);

/*
void VectorFromArrayFunction(Encapsulation *arr, u32 count, void *userdata)
{
   Real::FuncArray f = Real::FuncArray(userdata);
   vector<Encapsulation&> vec(vector<Encapsulation&>::iterator(arr),
                         vector<Encapsulation&>::iterator(arr + count));
   f(vec);
}

Real::Real(Real::FuncVector vecfunc, vector<Real&> vec)
{
   i32 count = vec.size();
   RealObject **parr = new RealObject* [count];
   for (int i=0;i<count;++i) {
      parr[i] = arr[i].m_pObject;
      parr[i]->AddRef();
   }
   RealArray *oarr = new RealArray(VectorFromArrayFunction, parr, count, (void*)(vecfunc));

   for (int i=0;i<count;++i)
      arr[i] = Real(new RealArrayElement(oarr, i);
   m_pObject = oarr;
}
*/
// operations. they create suitable RealUnary or RealBinary 
// object that link to the arguments

Real Real::operator - () const
{
   return Real(UnaryMinus, (const Real &) (*this), 0);
}

Real operator + (const Real &lhs, const Real &rhs)
{
    return Real(Plus, lhs, rhs, 0);
}

Real operator - (const Real &lhs, const Real &rhs)
{
    return Real(Minus, lhs, rhs, 0);
}

Real operator * (const Real &lhs, const Real &rhs)
{
    return Real(Multiply, lhs, rhs, 0);
}

Real operator / (const Real &lhs, const Real &rhs)
{
    return Real(Divide, lhs, rhs, 0);
}

Real operator * (const Real &lhs, int rhs)
{
    return Real(Multiply, lhs, UserInt(rhs));
}

Real operator / (const Real &lhs, int rhs)
{
    return Real(Divide, lhs, UserInt(rhs));
}

Real recip(const Real &rhs)
{
    return Real(recip, rhs, 0);
}

CreateRealConstant(Pi, pi)
CreateRealConstant(Ln2, ln2)
CreateUnaryRealFunction(sq)
CreateUnaryRealFunction(abs)
CreateUnaryRealFunction(sqrt)
CreateUnaryRealFunction(rsqrt)
CreateUnaryRealFunction(log)
CreateUnaryRealFunction(exp)
CreateBinaryRealFunction(atan2)
CreateUnaryRealFunction(tan)
CreateUnaryRealFunction(cos)
CreateUnaryRealFunction(sin)
CreateUnaryRealFunction(acos)
CreateUnaryRealFunction(asin)
CreateUnaryRealFunction(atan)
/*
Real rsqrt(const Real &rhs)
{
    return Real(rsqrt, rhs);
}

Real sqrt(const Real &rhs)
{
    return Real(sqrt, rhs);
}

Real log(const Real &rhs)
{
    return Real(log, rhs);
}

Real exp(const Real &rhs)
{
    return Real(exp, rhs);
}

Real abs(const Real &rhs)
{
    return Real(abs, rhs);
}

Real cos(const Real &rhs)
{
    return Real(cos, rhs);
}

Real sin(const Real &rhs)
{
    return Real(sin, rhs);
}

Real tan(const Real &rhs)
{
    return Real(tan, rhs);
}

Real acos(const Real &rhs)
{
    return Real(acos, rhs);
}

Real asin(const Real &rhs)
{
    return Real(asin, rhs);
}

Real atan(const Real &rhs)
{
    return Real(asin_primary, rhs);
}

Real atan2(const Real &l, const Real &r)
{
    return Real(atan2, l, r);
}

/*
Real rpi()
{
    return Real(Constants.rpi.m_pVal);
}


Real pi()
{
    return Real(pi);
}


Real ln2()
{
    return Real(Constants.ln2.m_pVal);
}

Real ln10()
{
    return Real(Constants.ln10.m_pVal);
}
*/

// the tricky part

struct RecursionSavedState {
    RealObject *ptr;
    int sibling;
};

Encapsulation PerformEvaluation(RealObject *pObj)
{
    assert(pObj);
    int maxpos = 0;
    if (!(pObj->HasEstimate() || pObj->m_Depth < EvaluationDepth)) {
        // depths range from 0 to this object's depth, inclusive
        RecursionSavedState* Stack = new RecursionSavedState[pObj->m_Depth+1];
        int sPos = 0;
        RealObject *pSib;
        pObj->AddRef();
        Stack[sPos].ptr = pObj;
        Stack[sPos].sibling = 0;
        
        try {
        // while we have not exhausted the depth-first search
        while (sPos >= 0) {
            // while there are more siblings of the currect object
            while (!!(pSib = Stack[sPos].ptr->GetSibling(Stack[sPos].sibling))) {
                ++Stack[sPos].sibling;
                // if we don't already have an Encapsulation
                if (!pSib->HasEstimate() /*&& pSib->m_Depth > EvaluationDepth*/) {
                    // save and reference the new child and move to it
                    ++sPos;
                    if (sPos > maxpos) maxpos = sPos;
                    pSib->AddRef();
                    Stack[sPos].ptr = pSib;
                    Stack[sPos].sibling = 0;
                }
            }
            // no more children of this node: prepare an Encapsulation
            // note we're now sure there will be no recursion here
            // because we've prepared estimates in all siblings
            Stack[sPos].ptr->GetEstimate();
            Stack[sPos].ptr->Release();
        
            // move up in the stack
            --sPos;
        }
        } catch (RealLibException e) {
            while (sPos>=0) {
                Stack[sPos--].ptr->Release();
            }
            delete [] Stack;
            throw;
        }

        delete [] Stack;
    } 
    
    // PerformEvaluation should not count as a use of the object
   // because the user working with Real's might request a value more than once
    ++pObj->m_EstimateRefs;
    return pObj->GetEstimate(); 
}

void PerformDeletion(RealObject *pObj)
{
    assert(pObj);
    if (pObj->GetRefCount() == 1 && pObj->m_Depth >= EvaluationDepth) {
        // depths range from 0 to this object's depth, inclusive
        RecursionSavedState* Stack = new RecursionSavedState[pObj->m_Depth+1];
        int sPos = 0;
        RealObject *pSib;
        Stack[sPos].ptr = pObj;
        Stack[sPos].sibling = 0;
        
        // while we have not exhausted the depth-first search
        while (sPos >= 0) {
            // while there are more siblings of the currect object
            while (!!(pSib = Stack[sPos].ptr->GetSibling(Stack[sPos].sibling))) {
                ++Stack[sPos].sibling;
                // if we have to delete it also
                if (pSib->GetRefCount() == 1 /*&& pSib->m_Depth > EvaluationDepth*/) {
                    // save the new child and move to it
                    // no need to reference, we're stealing the existing one
                    ++sPos;
                    Stack[sPos].ptr = pSib;
                    Stack[sPos].sibling = 0;
                } else pSib->Release();
            }
            // no more children of this node: delete it
            // do not release as this will trigger recursive
            // calls which we have already performed
            // note we will reach this point only for object that
            // really need to be deleted.
            assert(Stack[sPos].ptr->GetRefCount() == 1);
            Stack[sPos].ptr->NonRecursiveRelease();
        
            // move up in the stack
            --sPos;
        }
        
        delete [] Stack;
    } else
        pObj->Release();
}

/* obsolete older version

// evaluation can go very deep, and the stack is severely limited in depth
// to avoid problems, whenever Real evaluation or destruction should occur
// we first separate the tree to smaller segments, and then build up on
// them to perform the operation

// perform evaluation: build a list of all the objects that can be 
// reached in EvaluationDepth steps from root. Then build a list of
// all objects that can be reached from there, and so on till we reach
// the leaves. Then evaluate the lowest objects first -- they will remember
// the result so that evaluation on higher levels won't dig deeper than
// them.
Encapsulation PerformEvaluation(RealObject *pObj)
{
    #ifdef LIMIT_RECURSION
    // there's nothing to do if we have Encapsulation or if the depth is low
    if (!(pObj->HasEstimate() || pObj->m_Depth < EvaluationDepth)) {
        // start from the root
        EvalList *pFirst = new EvalList;
      pFirst->obj = pObj->AddRef();
        pFirst->next = 0;
        EvalList *pLast = pFirst;
        EvalList *pList, *pNewList;

        // grab the objects to pNewList
        pNewList = pLast;
        pLast->used = true;
        pLast->obj->GetDepthList(EvaluationDepth, FUNCTION(RealObject::HasEstimateReference), 0, &pNewList, true);
            // HasEstimateReference addresses two issues:
            //    if an object already has an Encapsulation, we don't need to go deeper
            //    if not, mark the object as visited so that we don't walk the
            //         full breadth of the tree (which is enormous)

        do {
            // attach to main list, reversed
            while (!pNewList->used) {
                pList = pNewList->next;
                pNewList->next = pFirst;
                pFirst = pNewList;
                pNewList = pList;
            }

            // work on main list
            while (!pFirst->used) {
                pLast = pLast->next = pFirst;
                pFirst = pFirst->next;
                pLast->next = 0;
                pLast->used = true;

                pLast->obj->m_EstimateRefs = 0;
                pLast->obj->GetDepthList(EvaluationDepth, FUNCTION(RealObject::HasEstimateReference), 0, &pNewList, true);
            }
        } while (!pNewList->used);

        // reverse list
        pLast = pFirst;
        pFirst = pFirst->next;
        pLast->next = 0;
        while (pFirst) {
            EvalList *pList = pFirst->next;
            pFirst->next = pLast;
            pLast = pFirst;
            pFirst = pList;
        }
        pFirst = pLast;

        // follow the list and evaluate from the bottom up
        while (pFirst) {
            EvalList *pList = pFirst;
            pFirst = pFirst->next;

            pList->obj->GetEstimate();   // this will not affect the number of queries to the Encapsulation 
            pList->obj->Release();           // because there is an extra reference, and release does not decrease EstimateRefs
            delete pList;
        }
    }
    #endif // LIMIT_RECURSION

    // PerformEvaluation should not count as a use of the object
   // because the user working with Real's might request a value more than once
    ++pObj->m_EstimateRefs;
    return pObj->GetEstimate();
}

// perform deletion: do it recursively only up to EvaluationDepth. When it is
// reached, remember the object and return. Then delete the remembered
// objects.
void PerformDeletion(RealObject *pRelObj)
{
    #ifdef LIMIT_RECURSION
    // if root has more references we don't need to delete
    if (pRelObj->GetRefCount() == 1 && pRelObj->m_Depth >= EvaluationDepth) {
        // start from root
        EvalList *pFirst = new EvalList;
      pFirst->obj = pRelObj;
        pFirst->next = 0;
        EvalList *pLast = pFirst;

        pLast->used = true;
        pLast->obj->GetDepthList(EvaluationDepth, FUNCTION(RealObject::StartRelease), FUNCTION(RealObject::FinishRelease), &pFirst, false);
            // StartRelease decreases the reference count and stops if
            // object wouldn't be deleted.
            // FinishRelease destroys the object if needed, but does not
            // release the siblings; they have already been released, or
            // have references in the EvalList and are scheduled to be.

        while (!pFirst->used) {
            pLast = pFirst;
            pFirst = pFirst->next;

            pLast->obj->GetDepthList(EvaluationDepth, FUNCTION(RealObject::StartRelease), FUNCTION(RealObject::FinishRelease), &pFirst, false);
            delete pLast;
        }
        delete pFirst;
    } else
    #endif // LIMIT_RECURSION
       pRelObj->Release();
}
*/

Real::~Real()
{
    PerformDeletion(m_pObject);
}

Real& Real::operator = (const Real &rhs)
{
    PerformDeletion(m_pObject);
    m_pObject = rhs.m_pObject;
    assert(m_pObject);
    m_pObject->AddRef();
    return *this;
}

#define PRECMULTIPLIER  2
#define PRECDIVISOR     1

bool IsProbableZero(const Encapsulation &r, u32 prec)
{
    // "probable zero": 0.0 is in the range of possible values,
    // and the current Encapsulation is smaller than the wanted precision
   return !r.IsNonZero() && r.weak_normalize() < -i32(prec);
}

#define DELETE(x) { if (x) { delete x; x = NULL; } }


// when resetting from an internal function, we must make sure
// Encapsulation's begin and finish functions are called correctly
// also exceptions from within MakeGoodEstimate/MakeNonZero
// have to call the finish function
unsigned InternalReset(unsigned precStart)
{
    Encapsulation::FinishComputation();
    unsigned precMax = g_precMax;
    unsigned res = FinalizeRealLib();
    InitializeRealLib(precStart, precMax, (g_memReq + precStart - 1) / precStart);
    Encapsulation::BeginComputation();
    return res;
}


// evaluate object.
// gradually increase the precision until the result can be considered
// correct to the requirements
Encapsulation MakeGoodEstimate(RealObject *pObj, i32 precision, i32 probzeroprec)
{
   if (g_WorkingPrecision*32 <= precision + 32) InternalReset(precision/32 + 3);
   u32 prec = g_WorkingPrecision;

   do {
      try {
         Encapsulation r(PerformEvaluation(pObj));
            if (r.GetRelativeError() < precision && !IsProbableZero(r, probzeroprec)) {
                throw PrecisionException("MakeGoodEstimate()");
            }
            return r;
      } catch (PrecisionException e) {
       if (prec >= g_precMax) throw;
            else InternalReset(prec = min(g_precMax, prec * PRECMULTIPLIER / PRECDIVISOR));
      } 
   } while (1);
}

Encapsulation MakeNonZero(RealObject *pObj)
{
   u32 prec = g_WorkingPrecision;

   do {
      try {
            Encapsulation r(PerformEvaluation(pObj));
            if (!r.IsNonZero()) {
                throw PrecisionException("MakeNonZero()");
            }
         return r;
      } catch (PrecisionException e) {
       if (prec >= g_precMax) throw;
            else InternalReset(prec = min(g_precMax, prec * PRECMULTIPLIER / PRECDIVISOR));
      }
   } while (1);
}

double Real::AsDouble() const
{
   Encapsulation::Computation comp;
    Encapsulation r(MakeGoodEstimate(m_pObject, 55, 1024));
    return r.weak_AsDouble();
}

char *Real::AsDecimal(char *buffer, unsigned lenwanted) const
{
   Encapsulation::Computation comp;
    u32 bitswanted = u32(ceil(LOG_2_10 * lenwanted));
    Encapsulation r(MakeGoodEstimate(m_pObject, bitswanted, bitswanted*2));

    if (IsProbableZero(r, bitswanted*2)) strcpy(buffer, "probable zero");
    else r.weak_AsDecimal(buffer, lenwanted);
   return buffer;
}

// looks like a strange function, and it is, because it
// can only return true, since equality is not decidable
// this function can also be exited via exceptions--
// which would mean that the number is zero up to the
// current max precision
bool Real::IsNonZero() const
{
   Encapsulation::Computation comp;
    MakeNonZero(m_pObject);
   return true;
}

bool Real::IsNegative() const
{
   Encapsulation::Computation comp;
    Encapsulation r(MakeNonZero(m_pObject));
    return r.IsNegative();
}

bool Real::IsPositive() const
{
   Encapsulation::Computation comp;
    Encapsulation r(MakeNonZero(m_pObject));
    return r.IsPositive();
}

std::ostream& operator << (std::ostream& out, const Real &re)
{
    Encapsulation::Computation comp;
    u32 pr = out.precision();
    if (out.flags() & std::ios::scientific) ++pr;
    u32 bitswanted = u32(ceil(LOG_2_10 * (pr)));
    u32 pzbits = bitswanted;
    if (!(out.flags() & std::ios::fixed)) pzbits *= 2;
    Encapsulation r(MakeGoodEstimate(re.m_pObject, bitswanted, pzbits));

    if (out.flags() & std::ios::fixed) {
        if (r.weak_ge(1)) {
            bitswanted += r.weak_normalize();
            r = MakeGoodEstimate(re.m_pObject, bitswanted, pzbits);
        }
    }

    if (IsProbableZero(r, pzbits) && !(out.flags() & std::ios::fixed)) out << "probable zero";
    else out << r;
   return out;
}

std::istream& operator >> (std::istream &in, Real &r)
{
    std::string s;
    in >> s;
    if (in) r = Real(s.c_str());
    
    return in;
}

} // namespace
