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

#include <stdlib.h>
#include <limits.h>
#include <exception>
#include "defs.h"
#include "LongFloat.h"
#include "RealEstimate.h"
#include "MachineEstimate.h"

namespace RealLib {

typedef ptrdiff_t UserInt;

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
    MachineEstimate m_mach;
    Estimate m_est;

public:
    Encapsulation(const Estimate &est)
        : m_est(est) {}
    Encapsulation(const MachineEstimate &mach)
        : m_mach(mach) {}

public:
    static void BeginComputation()
    {   if (UsingMachinePrecision)
        MachineEstimate::BeginComputation(); }  
    static void FinishComputation()
    {   if (UsingMachinePrecision)
        MachineEstimate::FinishComputation(); }

    class Computation {
    public:
        Computation()
        { BeginComputation(); }
        ~Computation()
        { FinishComputation(); }
    };

   Encapsulation(double v = 1.0)
   { if (UsingMachinePrecision)
        m_mach = MachineEstimate(v);
   else m_est = Estimate(v); }

   Encapsulation(const char *val)
   { if (UsingMachinePrecision)
        m_mach = MachineEstimate(val);
   else m_est = Estimate(val); }

    // error functions
   Encapsulation GetError() const
   { if (UsingMachinePrecision)
        return m_mach.GetError();
   else return m_est.GetError(); }

   Encapsulation& SetError(const Encapsulation &err)
   { if (UsingMachinePrecision)
        m_mach.SetError(err.m_mach);
   else m_est.SetError(err.m_est); 
   return *this; }

   Encapsulation& AddError(const Encapsulation &err)
   { if (UsingMachinePrecision)
        m_mach.AddError(err.m_mach);
   else m_est.AddError(err.m_est); 
   return *this; }
   

   // a lower bound on the correct binary digits
   // uses the exponents of the value and error to calculate it quickly
   i32 GetRelativeError() const 
   { if (UsingMachinePrecision)
        return m_mach.GetRelativeError();
   else return m_est.GetRelativeError(); }

   // get a rough estimate of the precision
   // used to determine the length of the approximations to functions
   u32 GetPrecision() const
   { if (UsingMachinePrecision)
        return m_mach.GetPrecision();
   else return m_est.GetPrecision(); }

   Encapsulation& SetPrecision(u32 prec)
;//      { m_Value.SetPrecision(prec); 
   //     return *this; }
      
   // comparisons
   // these come in two flavors, strong (true if real is in relation to rhs)
   bool IsPositive() const
   { if (UsingMachinePrecision)
        return m_mach.IsPositive();
   else return m_est.IsPositive(); }

   bool IsNegative() const
   { if (UsingMachinePrecision)
        return m_mach.IsNegative();
   else return m_est.IsNegative(); }

   bool IsNonZero() const
   { if (UsingMachinePrecision)
        return m_mach.IsNonZero();
   else return m_est.IsNonZero(); }


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
        return m_mach.weak_IsPositive();
   else return m_est.weak_IsPositive(); }

   bool weak_IsNegative() const
   { if (UsingMachinePrecision)
        return m_mach.weak_IsNegative();
   else return m_est.weak_IsNegative(); }
   
   bool weak_lt(const Encapsulation &rhs) const
   { if (UsingMachinePrecision)
        return m_mach.weak_lt(rhs.m_mach);
   else return m_est.weak_lt(rhs.m_est); }
   bool weak_eq(const Encapsulation &rhs) const
   { if (UsingMachinePrecision)
        return m_mach.weak_eq(rhs.m_mach);
   else return m_est.weak_eq(rhs.m_est); }

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
        return m_mach.weak_round();
   else return m_est.weak_round(); }

   // weak normalize, i.e. return an exponent such that 
   // a >> a.weak_normalize()
   // is in the range [0.5, 1).
   i32 weak_normalize() const
   { if (UsingMachinePrecision)
        return m_mach.weak_normalize();
   else return m_est.weak_normalize(); }
  
   // weak conversion
   double weak_AsDouble() const
   { if (UsingMachinePrecision)
        return m_mach.weak_AsDouble();
   else return m_est.weak_AsDouble(); }

   // output
   char *weak_AsDecimal(char *buffer, u32 buflen) const
   { if (UsingMachinePrecision)
        return m_mach.weak_AsDecimal(buffer, buflen);
   else return m_est.weak_AsDecimal(buffer, buflen); }
   
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
        return m_mach << howmuch;
    else return m_est << howmuch; }

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

/*
// shorthands
static inline 
Encapsulation operator * (i32 lhs, const Encapsulation &rhs)
{ return rhs * lhs; }
/*{ if (UsingMachinePrecision)
       return rhs.m_mach * lhs;
  else return rhs.m_est * lhs; }*/

/*
static inline 
Encapsulation operator / (i32 lhs, const Encapsulation &rhs)
{ return recip(rhs) * lhs; }
/*{ if (UsingMachinePrecision)
       return recip(rhs.m_mach) * lhs;
  else return recip(rhs.m_est) * lhs; }*/


// operations
static inline
Encapsulation UnaryMinus (const Encapsulation &arg, UserInt i)
{ if (UsingMachinePrecision)
       return -arg.m_mach;
  else return -arg.m_est; }

static inline
Encapsulation recip(const Encapsulation &arg, UserInt i)
{ if (UsingMachinePrecision)
       return recip(arg.m_mach);
  else return recip(arg.m_est); }

static inline
Encapsulation Plus (const Encapsulation &lhs, const Encapsulation &rhs, UserInt i)
{ if (UsingMachinePrecision)
       return lhs.m_mach + rhs.m_mach;
  else return lhs.m_est + rhs.m_est; }

static inline
Encapsulation Minus (const Encapsulation &lhs, const Encapsulation &rhs, UserInt i)
{ if (UsingMachinePrecision)
       return lhs.m_mach - rhs.m_mach;
  else return lhs.m_est - rhs.m_est; }

static inline
Encapsulation Multiply (const Encapsulation &lhs, const Encapsulation &rhs, UserInt i)
{ if (UsingMachinePrecision)
       return lhs.m_mach * rhs.m_mach;
  else return lhs.m_est * rhs.m_est; }

static inline
Encapsulation Divide (const Encapsulation &lhs, const Encapsulation &rhs, UserInt i)
{ if (UsingMachinePrecision)
       return lhs.m_mach / rhs.m_mach;
  else return lhs.m_est / rhs.m_est; }


// fast multiplication
static inline
Encapsulation Multiply (const Encapsulation &lhs, UserInt rhs)
{ if (UsingMachinePrecision)
       return lhs.m_mach * i32(rhs);
  else return lhs.m_est * i32(rhs); }

// and division
static inline
Encapsulation Divide (const Encapsulation &lhs, UserInt rhs)
{ if (UsingMachinePrecision)
       return lhs.m_mach / i32(rhs);
  else return lhs.m_est / i32(rhs); }



// C++-style output
static inline
std::ostream& operator <<(std::ostream &os, const Encapsulation &e)
{ if (UsingMachinePrecision)
       return os << e.m_mach;
  else return os << e.m_est; }

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
