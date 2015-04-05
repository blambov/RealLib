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

/*

  RealEncapsulation.h

  Encapsulates the different kinds of Estimates that can be used
  in the system so that RealObjects can refer to them in an uniform way.
  Classes:

    Encapsulation - an Estimate object with each method amounting to a
        switch on g_WorkingPrecision to choose a MachineEstimate or
        Estimate implementation
 */

#ifndef FILE_REAL_Encapsulation_H
#define FILE_REAL_Encapsulation_H

#define MachineEstimatePrecision 4
#define UsingMachinePrecision g_WorkingPrecision == MachineEstimatePrecision

#define REALLIB_ENCAPPTR_ALLOCS_MEMORY
// non-memory allocating version does not work with arrays!
// arrays should be redone to use the classes mechanism (almost as efficient, less error prone)

#include <stdlib.h>
#include <limits.h>
#include <assert.h>
#include <malloc.h>
#include <exception>
#include "defs.h"
#include "LongFloat.h"
#include "RealEstimate.h"
#include "MachineEstimate.h"

namespace RealLib {

typedef long long UserInt;

// class Encapsulation's definitions start here
class Encapsulation;

// operations
static inline   Encapsulation UnaryMinus (const Encapsulation &arg, UserInt user);
static inline   Encapsulation recip(const Encapsulation &arg, UserInt user);

static inline   Encapsulation Plus (const Encapsulation &lhs, const Encapsulation &rhs, UserInt user);
static inline   Encapsulation Minus (const Encapsulation &lhs, const Encapsulation &rhs, UserInt user);
static inline   Encapsulation Multiply (const Encapsulation &lhs, const Encapsulation &rhs, UserInt user);
static inline   Encapsulation Divide (const Encapsulation &lhs, const Encapsulation &rhs, UserInt user);

// fast multiplication
static inline   Encapsulation Multiply (const Encapsulation &lhs, UserInt rhs);
// and division
static inline   Encapsulation Divide (const Encapsulation &lhs, UserInt rhs);

static inline
std::ostream& operator <<(std::ostream &os, const Encapsulation &e);

class Encapsulation {
private:
public:
    // assuming MachineEstimate is as big as MachineEstimateBaseType
    MachineEstimateBaseType Storage[(sizeof (Estimate) + sizeof(MachineEstimateBaseType)-1)
                                    /sizeof(MachineEstimateBaseType)];

    Estimate& rwEstimate()
    { return *((Estimate*)(Storage)); }
    const Estimate& roEstimate() const
    { return *((const Estimate*)(Storage)); }

    MachineEstimate& rwMachineEstimate()
    { return *((MachineEstimate*)(Storage)); }
    const MachineEstimate& roMachineEstimate() const
    { return *((const MachineEstimate*)(Storage)); }

    //MachineEstimate m_mach;
    //Estimate m_est;

public:
    Encapsulation(const Estimate &est)
{ assert(!(UsingMachinePrecision));
new (Storage) Estimate(est); }
    Encapsulation(const MachineEstimate &mach)
    { assert(UsingMachinePrecision);
    new (Storage) MachineEstimate(mach); }

    Encapsulation(const Encapsulation &rhs)
    { if (UsingMachinePrecision)
        new (Storage) MachineEstimate(rhs.roMachineEstimate());
    else new (Storage) Estimate(rhs.roEstimate()); }

    Encapsulation& operator=(const Encapsulation &rhs)
    { if (UsingMachinePrecision)
        rwMachineEstimate() = rhs.roMachineEstimate();
    else rwEstimate() = rhs.roEstimate();  
    return *this; }

    ~Encapsulation()
    { if (UsingMachinePrecision)
        (rwMachineEstimate()).~MachineEstimate();
    else (rwEstimate()).~Estimate(); }

public:
    static void BeginComputation()
    {   if (UsingMachinePrecision)
        MachineEstimate::BeginComputation(); }  
    static void FinishComputation()
    {   if (UsingMachinePrecision)
        MachineEstimate::FinishComputation(); }

    bool IsValueValid() const
    { if (UsingMachinePrecision)
        return (roMachineEstimate()).IsValueValid();
    else return true; }

    class Computation {
    public:
        Computation()
    { BeginComputation(); }
        ~Computation()
        { FinishComputation(); }
    };

    Encapsulation(double v = 0.0)
    { if (UsingMachinePrecision)
        new (Storage) MachineEstimate(v);
    else new (Storage) Estimate(v); }

    Encapsulation(const char *val)
    { if (UsingMachinePrecision)
        new (Storage) MachineEstimate(val);
    else new (Storage) Estimate(val); }

    // error functions
    Encapsulation GetError() const
    { if (UsingMachinePrecision)
        return (roMachineEstimate()).GetError();
    else return (roEstimate()).GetError(); }

    Encapsulation& SetError(const Encapsulation &err)
    { if (UsingMachinePrecision)
        (rwMachineEstimate()).SetError(err.roMachineEstimate());
    else (rwEstimate()).SetError(err.roEstimate());
    return *this; }

    Encapsulation& AddError(const Encapsulation &err)
    { if (UsingMachinePrecision)
        (rwMachineEstimate()).AddError(err.roMachineEstimate());
    else (rwEstimate()).AddError(err.roEstimate());
    return *this; }


    // a lower bound on the correct binary digits
    // uses the exponents of the value and error to calculate it quickly
    i32 GetRelativeError() const
    { if (UsingMachinePrecision)
        return (roMachineEstimate()).GetRelativeError();
    else return (roEstimate()).GetRelativeError(); }

    // get a rough estimate of the precision
    // used to determine the length of the approximations to functions
    u32 GetPrecision() const
    { if (UsingMachinePrecision)
        return (roMachineEstimate()).GetPrecision();
    else return (roEstimate()).GetPrecision(); }

    Encapsulation& SetPrecision(u32 prec)
    ;//      { m_Value.SetPrecision(prec);
    //     return *this; }

    // comparisons
    // these come in two flavors, strong (true if real is in relation to rhs)
    bool IsPositive() const
    { if (UsingMachinePrecision)
        return (roMachineEstimate()).IsPositive();
    else return (roEstimate()).IsPositive(); }

    bool IsNegative() const
    { if (UsingMachinePrecision)
        return (roMachineEstimate()).IsNegative();
    else return (roEstimate()).IsNegative(); }

    bool IsNonZero() const
    { if (UsingMachinePrecision)
        return (roMachineEstimate()).IsNonZero();
    else return (roEstimate()).IsNonZero(); }


    // left like this... maybe should be changed to use appropriate versions.
    /*
    bool operator < (const Encapsulation &rhs) const
      { return (*this - rhs).IsNegative(); }   
   bool operator > (const Encapsulation &rhs) const
      { return (*this - rhs).IsPositive(); }   
   bool operator != (const Encapsulation &rhs) const
      { return (*this - rhs).IsNonZero(); }   
     */

    bool weak_IsPositive() const
    { if (UsingMachinePrecision)
        return (roMachineEstimate()).weak_IsPositive();
    else return (roEstimate()).weak_IsPositive(); }

    bool weak_IsNegative() const
    { if (UsingMachinePrecision)
        return (roMachineEstimate()).weak_IsNegative();
    else return (roEstimate()).weak_IsNegative(); }

    bool weak_lt(const Encapsulation &rhs) const
    { if (UsingMachinePrecision)
        return (roMachineEstimate()).weak_lt(rhs.roMachineEstimate());
    else return (roEstimate()).weak_lt(rhs.roEstimate()); }
    bool weak_eq(const Encapsulation &rhs) const
    { if (UsingMachinePrecision)
        return (roMachineEstimate()).weak_eq(rhs.roMachineEstimate());
    else return (roEstimate()).weak_eq(rhs.roEstimate()); }

    bool weak_gt(const Encapsulation &rhs) const
    { return rhs.weak_lt(*this); }

    bool weak_le(const Encapsulation &rhs) const
    { return !weak_gt(rhs); }
    bool weak_ne(const Encapsulation &rhs) const
    { return !weak_eq(rhs); }
    bool weak_ge(const Encapsulation &rhs) const
    { return !weak_lt(rhs); }

    // among the weak operations is also rounding
    // the returned Encapsulation is assumed exact
    // only to be used on periodic functions!
    Encapsulation weak_round() const
    { if (UsingMachinePrecision)
        return (roMachineEstimate()).weak_round();
    else return (roEstimate()).weak_round(); }

    // weak normalize, i.e. return an exponent such that
    // a >> a.weak_normalize()
    // is in the range [0.5, 1).
    i32 weak_normalize() const
    { if (UsingMachinePrecision)
        return (roMachineEstimate()).weak_normalize();
    else return (roEstimate()).weak_normalize(); }

    // weak conversion
    double weak_AsDouble() const
    { if (UsingMachinePrecision)
        return (roMachineEstimate()).weak_AsDouble();
    else return (roEstimate()).weak_AsDouble(); }

    // output
    char *weak_AsDecimal(char *buffer, u32 buflen) const
    { if (UsingMachinePrecision)
        return (roMachineEstimate()).weak_AsDecimal(buffer, buflen);
    else return (roEstimate()).weak_AsDecimal(buffer, buflen); }

    friend   Encapsulation UnaryMinus (const Encapsulation &arg, UserInt user);
    friend   Encapsulation recip(const Encapsulation &arg, UserInt user);

    friend   Encapsulation Plus (const Encapsulation &lhs, const Encapsulation &rhs, UserInt user);
    friend   Encapsulation Minus (const Encapsulation &lhs, const Encapsulation &rhs, UserInt user);
    friend   Encapsulation Multiply (const Encapsulation &lhs, const Encapsulation &rhs, UserInt user);
    friend   Encapsulation Divide (const Encapsulation &lhs, const Encapsulation &rhs, UserInt user);

    // fast multiplication
    friend   Encapsulation Multiply (const Encapsulation &lhs, UserInt rhs);
    // and division
    friend   Encapsulation Divide (const Encapsulation &lhs, UserInt rhs);

    // binary shift
    Encapsulation operator << (i32 howmuch) const
    { if (UsingMachinePrecision)
        return (roMachineEstimate()) << howmuch;
    else return (roEstimate()) << howmuch; }

    Encapsulation operator >> (i32 howmuch) const
    { return *this << -howmuch; }

    /*
    Encapsulation& operator += (const Encapsulation &rhs)
       { return *this = *this + rhs; }
    Encapsulation& operator -= (const Encapsulation &rhs)
       { return *this = *this - rhs; }
    Encapsulation& operator *= (const Encapsulation &rhs)
       { return *this = *this * rhs; }
    Encapsulation& operator /= (const Encapsulation &rhs)
       { return *this = *this / rhs; }


    Encapsulation& operator >>= (i32 rhs) 
       { return *this = *this >> rhs; }
    Encapsulation& operator <<= (i32 rhs) 
       { return *this = *this << rhs; }
    Encapsulation& operator *= (i32 rhs) 
       { return *this = *this * rhs; }
    Encapsulation& operator /= (i32 rhs) 
       { return *this = *this / rhs; }*/

    // should probably be somewhere else
    // conversion to string
    // char *AsDecimal(char *buffer, u32 buflen);
    friend  
    std::ostream& operator <<(std::ostream &os, const Encapsulation &e);

};

class EncapsulationPointer {
private:
    Encapsulation *ptr;
#ifndef REALLIB_ENCAPPTR_ALLOCS_MEMORY
#ifndef REALLIB_MACHINEESTIMATE_ALIGNMENT_REQUIRED
    MachineEstimateBaseType Storage[(sizeof (Estimate) + sizeof(MachineEstimateBaseType)-1)
                                    /sizeof(MachineEstimateBaseType)];
#else
    MachineEstimateBaseType Storage[(sizeof (Estimate) + sizeof(MachineEstimateBaseType)-1)
                                    /sizeof(MachineEstimateBaseType)+1];
#endif
#endif

    explicit EncapsulationPointer(Encapsulation *p)
    : ptr(p) { assert(p); };
public:
    static EncapsulationPointer
    FromPointer(Encapsulation *p)
    { return EncapsulationPointer(p); }


    // copy protection
    //EncapsulationPointer(const EncapsulationPointer &rhs);

public:
    EncapsulationPointer()
: ptr(NULL) {}

    EncapsulationPointer(const Encapsulation &val)

#ifdef REALLIB_ENCAPPTR_ALLOCS_MEMORY
#ifdef REALLIB_MACHINEESTIMATE_ALIGNMENT_REQUIRED
    : ptr((Encapsulation*)_mm_malloc(sizeof(Encapsulation), sizeof(MachineEstimateBaseType)))
    {
        if (!ptr) throw std::bad_alloc();
        else new (ptr) Encapsulation(val);
    }
#else
    : ptr(new Encapsulation(val))
          {}
#endif
#else
#ifdef REALLIB_MACHINEESTIMATE_ALIGNMENT_REQUIRED
    : ptr((Encapsulation*)(((ptrdiff_t)Storage)+15 & ~15))
#else
            : ptr((Encapsulation*)Storage)
#endif
              {
        new(ptr) Encapsulation(val);
              }
#endif

    EncapsulationPointer(unsigned ArraySize)
#ifdef REALLIB_MACHINEESTIMATE_ALIGNMENT_REQUIRED
    : ptr((Encapsulation*)_mm_malloc(sizeof(Encapsulation)*ArraySize, sizeof(MachineEstimateBaseType)))
    {
        if (!ptr) throw std::bad_alloc();
        else for (int i=0;i<ArraySize;++i) {
            new (ptr+i) Encapsulation;
        }
    }
#else
    : ptr(new Encapsulation[ArraySize])
          {}
#endif

    ~EncapsulationPointer()
    {}

    EncapsulationPointer& operator=(const EncapsulationPointer &rhs)
#ifdef REALLIB_ENCAPPTR_ALLOCS_MEMORY
    { ptr = rhs.ptr; return *this; }
#else
    {
        //*ptr = *rhs.ptr; return *this; }
        if (!ptr)
#ifdef REALLIB_MACHINEESTIMATE_ALIGNMENT_REQUIRED
            ptr = ((Encapsulation*)(((ptrdiff_t)Storage)+15 & ~15));
#else
        ptr = ((Encapsulation*)Storage);
#endif

        if (rhs.ptr)
            *ptr = *rhs.ptr;
        else ptr = NULL;

        return *this; }
#endif

    void Release()
    {
#ifdef REALLIB_ENCAPPTR_ALLOCS_MEMORY
#ifdef REALLIB_MACHINEESTIMATE_ALIGNMENT_REQUIRED
        assert(ptr);
        ptr->~Encapsulation();
        _mm_free(ptr);
#else
        delete ptr;
#endif
#else
        ptr->~Encapsulation();
#endif
        ptr = NULL;
    }

    void ReleaseArray(unsigned ArraySize)
    {
#ifdef REALLIB_MACHINEESTIMATE_ALIGNMENT_REQUIRED
        assert(ptr);
        for (int i=0;i<ArraySize;++i)
            ptr[i].~Encapsulation();
        _mm_free(ptr);
        ptr = NULL;
#else
        delete [] ptr;
        ptr = NULL;
#endif
    }

    operator bool()
          { return !!ptr; }

    Encapsulation& operator *() { return *ptr; }
    const Encapsulation& operator *() const { return *ptr; }

    typedef Encapsulation *EncPtr;
    operator EncPtr() { return ptr; }

    EncapsulationPointer operator+ (i32 v)
    { return EncapsulationPointer(ptr+v); }
};


template <class DEST, class SRC>
static inline DEST explicit_cast(SRC arg)
{ return DEST(arg); }

#ifdef _MSC_VER
template <>
static inline 
Estimate explicit_cast<Estimate>(const Encapsulation &arg)
{ return arg.roEstimate(); }
#else
template <>
inline
Estimate explicit_cast<Estimate>(Encapsulation arg)
{ return arg.roEstimate(); }
#endif

#ifdef _MSC_VER
template <>
static inline 
MachineEstimate explicit_cast<MachineEstimate>(const Encapsulation &arg)
{ return arg.roMachineEstimate(); }
#else
template <>
inline
MachineEstimate explicit_cast<MachineEstimate>(Encapsulation arg)
{ return arg.roMachineEstimate(); }
#endif

/*
  template <>
  Estimate& ExplicitConvert()
  { return (roEstimate()); }
  template <>
  MachineEstimate& ExplicitConvert()
  { return (roMachineEstimate()); }*/


/*
// shorthands
static inline 
Encapsulation operator * (i32 lhs, const Encapsulation &rhs)
{ return rhs * lhs; }
/*{ if (UsingMachinePrecision)
       return rhs.roMachineEstimate() * lhs;
  else return rhs.roEstimate() * lhs; }*/

/*
static inline 
Encapsulation operator / (i32 lhs, const Encapsulation &rhs)
{ return recip(rhs) * lhs; }
/*{ if (UsingMachinePrecision)
       return recip(rhs.roMachineEstimate()) * lhs;
  else return recip(rhs.roEstimate()) * lhs; }*/


// operations
static inline
Encapsulation UnaryMinus (const Encapsulation &arg, UserInt i)
{ if (UsingMachinePrecision)
    return -arg.roMachineEstimate();
else return -arg.roEstimate(); }

static inline
Encapsulation recip(const Encapsulation &arg, UserInt i)
{ if (UsingMachinePrecision)
    return recip(arg.roMachineEstimate());
else return recip(arg.roEstimate()); }

static inline
Encapsulation Plus (const Encapsulation &lhs, const Encapsulation &rhs, UserInt i)
{ if (UsingMachinePrecision)
    return lhs.roMachineEstimate() + rhs.roMachineEstimate();
else return lhs.roEstimate() + rhs.roEstimate(); }

static inline
Encapsulation Minus (const Encapsulation &lhs, const Encapsulation &rhs, UserInt i)
{ if (UsingMachinePrecision)
    return lhs.roMachineEstimate() - rhs.roMachineEstimate();
else return lhs.roEstimate() - rhs.roEstimate(); }

static inline
Encapsulation Multiply (const Encapsulation &lhs, const Encapsulation &rhs, UserInt i)
{ if (UsingMachinePrecision)
    return lhs.roMachineEstimate() * rhs.roMachineEstimate();
else return lhs.roEstimate() * rhs.roEstimate(); }

static inline
Encapsulation Divide (const Encapsulation &lhs, const Encapsulation &rhs, UserInt i)
{ if (UsingMachinePrecision)
    return lhs.roMachineEstimate() / rhs.roMachineEstimate();
else return lhs.roEstimate() / rhs.roEstimate(); }


// fast multiplication
static inline
Encapsulation Multiply (const Encapsulation &lhs, UserInt rhs)
{ if (UsingMachinePrecision)
    return lhs.roMachineEstimate() * i32(rhs);
else return lhs.roEstimate() * i32(rhs); }

// and division
static inline
Encapsulation Divide (const Encapsulation &lhs, UserInt rhs)
{ if (UsingMachinePrecision)
    return lhs.roMachineEstimate() / i32(rhs);
else return lhs.roEstimate() / i32(rhs); }



// C++-style output
static inline
std::ostream& operator <<(std::ostream &os, const Encapsulation &e)
{ if (UsingMachinePrecision)
    return os << e.roMachineEstimate();
else return os << e.roEstimate(); }

// array interface

template <class TYPE, long offset = sizeof (TYPE)>
// GCC does not seem to support the above
//template <class TYPE, long offset>
class ArrayInterface {
    char *arr;
    long count;

public:
    ArrayInterface(TYPE *p, int c)
    : arr((char*)p), count(c) {}

    long size() { return count; }
    TYPE& operator[] (int index) 
    { return *((TYPE*)(arr + index*offset)); }
};

}   // namespace

#endif
