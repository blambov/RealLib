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
                Stack[sPos].ptr->Release(0);
                // Release(0) does not decrease m_EstimateRefs, but GetEstimate
                // does. We've added an extra reference, which is now taken
                // care of.

                // move up in the stack
                --sPos;
            }
        } catch (...) {
            // clear the recursion stack in case of exception
            while (sPos>=0) {
                Stack[sPos--].ptr->Release();
            }
            delete [] Stack;
            throw;
        }

        delete [] Stack;
    } 

    // PerformEvaluation should not count as a use of the object
    // because the user working with Reals might request a value more than once
    pObj->AddRef();
    Encapsulation res(pObj->GetEstimate()); 
    pObj->Release(0);
    return res;
    // Release(0) does not decrease the cached value reference count
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
        pObj->Release(1);
}

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
    return !r.IsNonZero() && Minus(absEncapsulation(r, 0) << prec, Encapsulation(1.0), 0).IsNegative(); 

    //    abs(r) < (Encapsulation(1) >> prec);
    //    && (r+r.GetError()).weak_normalize() < -i32(prec)
    //    && (r-r.GetError()).weak_normalize() < -i32(prec);
}

#define DELETE(x) { if (x) { delete x; x = NULL; } }

// when resetting from an internal function, we must make sure
// Encapsulation's begin and finish functions are called correctly
// also exceptions from within MakeGoodEstimate/MakeNonZero
// have to call the finish function
unsigned InternalReset(unsigned precStart)
{
    //cout << "Internal reset, prec " << g_WorkingPrecision << " new " << precStart << endl;
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
    //cout << "MakeGoodEstimate " << pObj << " prec " << precision << " pzprec " << probzeroprec << endl;
    if (g_WorkingPrecision*32 <= precision + 32) InternalReset(precision/32 + 3);
    u32 prec = g_WorkingPrecision;

    do {
        try {
            Encapsulation r(PerformEvaluation(pObj));
            if (!r.IsValueValid() || r.GetRelativeError() < precision && !IsProbableZero(r, probzeroprec)) {
                throw PrecisionException("MakeGoodEstimate()");
            }
            return r;
        } catch (PrecisionException &e) {
#ifdef REALLIB_SHOW_PRECISION_INCREASES
            using namespace std;
            cout << e;
            if (prec >= g_precMax) {
                cout << ", giving up.\n";
                throw;
            } else {
                cout << ", raising precision from " << prec;
                InternalReset(prec = min(g_precMax, prec * PRECMULTIPLIER / PRECDIVISOR));
                cout << " to " << prec << endl;
            }
#else
            e;
            if (prec >= g_precMax) throw;
            else InternalReset(prec = min(g_precMax, prec * PRECMULTIPLIER / PRECDIVISOR));
#endif
        }
    } while (1);
}

Encapsulation MakeNonZero(RealObject *pObj)
{
    u32 prec = g_WorkingPrecision;

    do {
        try {
            Encapsulation r(PerformEvaluation(pObj));
            if (!r.IsValueValid() || !r.IsNonZero()) {
                throw PrecisionException("MakeNonZero()");
            }
            return r;
        } catch (PrecisionException e) {
#ifdef REALLIB_SHOW_PRECISION_INCREASES
            using namespace std;

            cout << e;
            if (prec >= g_precMax) {
                cout << ", giving up.\n";
                throw;
            } else {
                cout << ", raising precision from " << prec;
                InternalReset(prec = min(g_precMax, prec * PRECMULTIPLIER / PRECDIVISOR));
                cout << " to " << prec << endl;
            }
#else
            if (prec >= g_precMax) throw;
            else InternalReset(prec = min(g_precMax, prec * PRECMULTIPLIER / PRECDIVISOR));
#endif
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

    if (IsProbableZero(r, bitswanted*2)) 
#ifdef _MSC_VER
#pragma warning (push)
#pragma warning (disable:4996)
        strcpy(buffer, "probable zero");
#pragma warning (pop)
#else
    strcpy(buffer, "probable zero");
#endif
    else r.weak_AsDecimal(buffer, lenwanted);
    return buffer;
}

// looks like a strange function, and it is, because it
// can only return true, since equality is not decidable.
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

    // important: r should only be alive while there's no chance of
    // reiteration
    {
        Encapsulation r(MakeGoodEstimate(re.m_pObject, bitswanted, pzbits));
        if (!(out.flags() & std::ios::fixed) || r.weak_lt(1)) {
            if (IsProbableZero(r, pzbits) && !(out.flags() & std::ios::fixed)) out << "probable zero";
            else out << r;
            return out;
        }
        bitswanted += r.weak_normalize();
    } 
    // r is destroyed here before the possible reiteration
    // which would not occur if r was good enough (due to caching)

    return out << MakeGoodEstimate(re.m_pObject, bitswanted, pzbits);
}

std::istream& operator >> (std::istream &in, Real &r)
{
    std::string s;
    in >> s;
    if (in) r = Real(s.c_str());

    return in;
}

} // namespace
