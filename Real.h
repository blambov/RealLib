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

  Real.h

  This file defines the class Real, which is the frontend of the
  library. It can be used like a normal built-in scalar, and
  provides useful conversions and tests. 

*/

#ifndef FILE_REAL_H
#define FILE_REAL_H

#include <vector>
#include <valarray>
#include "RealEncapsulation.h"
#include "RealFuncs.h"

namespace RealLib {

// initialize library. starts with precision precStart
// any operation that requires crossing the precMax boundary
// will abort with an exception
void InitializeRealLib(unsigned precStart = MachineEstimatePrecision, unsigned precMax = 100000, unsigned numEstAtStart = 1000);

// finalize. releases all memory in data objects. 
// returns the precision reached
unsigned FinalizeRealLib();

// reset calculations. useful if two separate calculations are to
// be carried out and the second wouldn't need as high precision.
// equivalent to finalize followed by initialize
unsigned ResetRealLib(unsigned precStart);

// returns the current precision
static inline unsigned GetCurrentPrecision()
{ return g_WorkingPrecision; }

// oracle functions take this form
typedef const char* (*OracleFunction) (unsigned precision);

// hide implementation details
class RealObject;
class Encapsulation;

class Real;

// operators
Real operator + (const Real &lhs, const Real &rhs);
Real operator - (const Real &lhs, const Real &rhs);
Real operator * (const Real &lhs, const Real &rhs);
Real operator / (const Real &lhs, const Real &rhs);

Real recip(const Real &arg);

Real operator * (const Real &lhs, int rhs);
Real operator / (const Real &lhs, int rhs);

// C++-style input/output
std::ostream& operator <<(std::ostream &out, const Real& r);
std::istream& operator >>(std::istream &in, Real &r);

// some definitions for array functions

// an array function must be declared like this:
// template <class TYPE, class ARRAY>
// void ArrayFunc(ARRAY &arr, void *otherData);
// 
// it must treat arr to have this interface:
// template <class TYPE>
// class ArrayTemplate {
// public:
//  long size();
//  TYPE& operator[] (long index);  // no bounds checking!
// };
typedef void (*RealFuncArray) (ArrayInterface<Encapsulation> &arr, UserInt userData);

template <class ARRAY>
void RealApplyArrayFunction(RealFuncArray arfunc, ARRAY &arr, UserInt user);

// the main object. can be used like a built-in type
class Real {
private:
    RealObject *m_pObject;

    // constructor
    Real(RealObject *pObject);

public:
    typedef Encapsulation (*FuncNullary) (unsigned int prec, UserInt userData); // real constant
    typedef Encapsulation (*FuncUnary) (const Encapsulation &arg, UserInt userData);
    typedef Encapsulation (*FuncBinary) (const Encapsulation &left, const Encapsulation &right, UserInt userData);
    //typedef Encapsulation (*FuncUnary) (const Encapsulation &arg);
    //typedef Encapsulation (*FuncBinary) (const Encapsulation &left, const Encapsulation &right);
    //typedef Encapsulation (*FuncBinaryOnInt) (const Encapsulation &left, long right);
    //typedef void (*FuncArray) (Encapsulation *array, unsigned int count, void *userData);
//   typedef Encapsulation (*FuncVector) (std::vector<Encapsulation&> vec);
   
   // copy constructor
    Real(const Real &src);

    // conversion from other types
    Real(const double src = 0);
    Real(const char *src);

   // reals created by operations
    Real(OracleFunction oracle);  // real constant by oracle, i.e. function returning char*
   Real(FuncNullary constant, UserInt user);   // real constant by nullary function. i.e. returning Encapsulation
   //Real(FuncUnary unfunc, const Real &rhs);
   //Real(FuncBinary binfunc, const Real &lhs, const Real &rhs);
   Real(FuncUnary unfunc, const Real &rhs, UserInt user);
   Real(FuncBinary binfunc, const Real &lhs, const Real &rhs, UserInt user);
   //Real(FuncBinaryOnInt binfunc, const Real &lhs, long rhs);
   
    template <class ARRAY>
    friend void RealApplyArrayFunction(RealFuncArray arfunc, ARRAY &arr, UserInt user);


    // assignment
    Real& operator = (const Real &rhs);

    // destructor
    ~Real();

    // conversion to other types
    double AsDouble() const;
    char* AsDecimal(char *buffer, unsigned lenwanted) const;

    // operators
    Real operator - () const;

    // comparison operator, semi-decidable
   bool IsNegative() const;
   bool IsPositive() const;
   bool IsNonZero() const;

    const Real& ForceNonZero() const
    { IsNonZero();
      return *this; }

    // shorthands
    Real& operator += (const Real &rhs)
    { return *this = *this + rhs; }
    Real& operator -= (const Real &rhs)
    { return *this = *this - rhs; }
    Real& operator *= (const Real &rhs)
    { return *this = *this * rhs; }
    Real& operator /= (const Real &rhs)
    { return *this = *this / rhs; }
    
    Real& operator *= (int rhs)
    { return *this = *this * rhs; }
    Real& operator /= (int rhs)
    { return *this = *this / rhs; }
    
    Real& operator *= (double rhs)
    { return *this = *this * Real(rhs); }
    Real& operator /= (double rhs)
    { return *this = *this / Real(rhs); }
    
    friend std::ostream& operator <<(std::ostream &out, const Real& r);
    friend std::istream& operator >>(std::istream &in, Real &r);

};

// to be able to work on the faster levels (MachineEstimate and Estimate),
// functions on Reals have to be declared as one of

/*
    // nullary real function (constant)
    template <class TYPE>
    TYPE name(unsigned int prec);

    // nullary real function with integer argument
    template <class TYPE>
    TYPE name(unsigned int prec, UserInt uint);

    // unary real function
    template <class TYPE>
    TYPE name(const TYPE &arg);

    // unary real function with integer argument
    template <class TYPE>
    TYPE name(const TYPE& arg, UserInt uint);

    // binary real function
    template <class TYPE>
    TYPE name(const TYPE &lhs, const TYPE &rhs);

    // binary real function with integer argument
    template <class TYPE>
    TYPE name(const TYPE& lhs, const TYPE &rhs, UserInt uint);

    // real function on array
    template <class TYPE, class ARRAY>
    TYPE name(ARRAY &arg);

    // real function on array with integer argument
    template <class TYPE>
    TYPE name(ARRAY &arg, UserInt int);

    // where Array has the following interface
    template <class TYPE>
    class ArrayInterface {
    public:
        long size();
        TYPE& operator[] (long index);  
    };
*/


#define CreateRealConstant(const_name, func_name) \
    Encapsulation func_name ## Encapsulation(unsigned int prec, UserInt user) { \
        if (UsingMachinePrecision) return func_name<MachineEstimate>(prec); \
        else return func_name<Estimate>(prec); \
    } \
    const Real const_name(func_name ## Encapsulation, 0);

#define CreateNullaryRealFunction(name) \
    Encapsulation name ## Encapsulation(unsigned int prec, UserInt user) { \
        if (UsingMachinePrecision) return name<MachineEstimate>(prec); \
        else return name<Estimate>(prec); \
    } \
    Real name() { \
        return Real(name ## Encapsulation, 0); \
    }

#define CreateIntRealFunction(name) \
    Encapsulation name ## Encapsulation(unsigned int prec, UserInt user) { \
        if (UsingMachinePrecision) return name<MachineEstimate>(prec, user); \
        else return name<Estimate>(prec, user); \
    } \
    Real name(UserInt user) { \
        return Real(name ## Encapsulation, user); \
    }

#define CreateUnaryRealFunction(name) \
    Encapsulation name ## Encapsulation(const Encapsulation &arg, UserInt user) { \
        if (UsingMachinePrecision) return name(arg.m_mach); \
        else return name(arg.m_est); \
    } \
    Real name(const Real &arg) { \
        return Real(name ## Encapsulation, arg, 0); \
    }

#define CreateUnaryAndIntRealFunction(name) \
    Encapsulation name ## Encapsulation(const Encapsulation &arg, UserInt user) { \
        if (UsingMachinePrecision) return name(arg.m_mach, user); \
        else return name(arg.m_est, user); \
    } \
    Real name(const Real &arg, UserInt rhs) { \
        return Real(name ## Encapsulation, arg, rhs); \
    }

#define CreateBinaryRealFunction(name) \
    Encapsulation name ## Encapsulation(const Encapsulation &lhs, const Encapsulation &rhs, UserInt user) { \
        if (UsingMachinePrecision) return name(lhs.m_mach, rhs.m_mach); \
        else return name(lhs.m_est, rhs.m_est); \
    } \
    Real name(const Real &arg, const Real &rhs) { \
        return Real(name ## Encapsulation, arg, rhs, 0); \
    }

#define CreateBinaryAndIntRealFunction(name) \
    Encapsulation name ## Encapsulation(const Encapsulation &lhs, const Encapsulation &rhs, UserInt user) { \
        if (UsingMachinePrecision) return name(lhs.m_mach, rhs.m_mach, user); \
        else return name(lhs.m_est, rhs.m_est, user); \
    } \
    Real name(const Real &arg, const Real &rhs, UserInt user) { \
        return Real(name ## Encapsulation, arg, rhs, user); \
    }

/*
#define CreateArrayAndOtherDataRealFunction(name) \
    void name ## Encapsulation(ArrayInterface<Encapsulation> &arr, void *otherData) { \
        if (UsingMachinePrecision) name(ArrayInterface<MachineEstimate, sizeof Encapsulation>(&arr[0].m_mach, arr.size()), otherData); \
        else name(ArrayInterface<Estimate, sizeof Encapsulation>(&arr[0].m_est, arr.size()), otherData); \
    } \
    Real name(Real *ptr, long count, void *otherData) { \
        return Real(name ## Encapsulation, ArrayInterface<Real>(ptr, count), otherData); \
    } \
    Real name(std::valarray<Real> &arr, long otherData) { \
        return Real(name ## Encapsulation, arr, otherData)); \
    } \
    Real name(std::vector<Real> &arr, long otherData) { \
        return Real(name ## Encapsulation, arr, otherData)); \
    }*/

#define CreateArrayAndIntRealFunction(name) \
    void name ## Encapsulation(ArrayInterface<Encapsulation> &arr, UserInt user) { \
        if (UsingMachinePrecision) { \
            ArrayInterface<MachineEstimate, sizeof (Encapsulation)> am(&arr[0].m_mach, arr.size()); \
            name<MachineEstimate, ArrayInterface<MachineEstimate, sizeof (Encapsulation)> >(am, user); \
        } else { \
            ArrayInterface<Estimate, sizeof (Encapsulation)> ae(&arr[0].m_est, arr.size()); \
            name<Estimate, ArrayInterface<Estimate, sizeof (Encapsulation)> > (ae, user);\
        } \
    } \
    void name(Real *ptr, long count, UserInt otherData) { \
        ArrayInterface<Real> arr(ptr, count); \
        RealApplyArrayFunction(name ## Encapsulation, arr, otherData); \
    } \
    void name(std::valarray<Real> &arr, UserInt otherData) { \
        RealApplyArrayFunction(name ## Encapsulation, arr, otherData); \
    } \
    void name(std::vector<Real> &arr, UserInt otherData) { \
        RealApplyArrayFunction(name ## Encapsulation, arr, otherData); \
    } 

#define CreateArrayRealFunction(name) \
    void name ## Encapsulation(ArrayInterface<Encapsulation> &arr, UserInt otherData) { \
        if (UsingMachinePrecision) { \
            ArrayInterface<MachineEstimate, sizeof (Encapsulation)> am(&arr[0].m_mach, arr.size()); \
            name<MachineEstimate, ArrayInterface<MachineEstimate, sizeof (Encapsulation)> >(am); \
        } else { \
            ArrayInterface<Estimate, sizeof (Encapsulation)> ae(&arr[0].m_est, arr.size()); \
            name<Estimate, ArrayInterface<Estimate, sizeof (Encapsulation)> > (ae);\
        } \
    } \
    void name(Real *ptr, long count) { \
        ArrayInterface<Real> arr(ptr, count); \
        RealApplyArrayFunction(name ## Encapsulation, arr, NULL); \
    } \
    void name(std::valarray<Real> &arr) { \
        RealApplyArrayFunction(name ## Encapsulation, arr, NULL)); \
    } \
    void name(std::vector<Real> &arr) { \
        RealApplyArrayFunction(name ## Encapsulation, arr, NULL)); \
    }


// constants
extern const Real Pi;
extern const Real Ln2;

// functions
Real abs(const Real &arg);
Real sqrt(const Real &arg);
Real rsqrt(const Real &arg);
Real log(const Real &arg);
Real exp(const Real &arg);
Real abs(const Real &arg);
Real atan2(const Real &y, const Real &x);
Real tan(const Real &arg);
Real cos(const Real &arg);
Real sin(const Real &arg);
Real acos(const Real &arg);
Real asin(const Real &arg);
Real atan(const Real &arg);

static inline Real operator * (int lhs, const Real &rhs)
{ return rhs * lhs; }
static inline Real operator / (int lhs, const Real &rhs)
{ return recip(rhs) * lhs; };

static inline Real operator * (const Real &lhs, double rhs)
{ return lhs * Real(rhs); }
static inline Real operator / (const Real &lhs, double rhs)
{ return lhs * Real(rhs); }
static inline Real operator * (double lhs, const Real &rhs)
{ return Real(lhs) * rhs; }
static inline Real operator / (double lhs, const Real &rhs)
{ return Real(lhs) / rhs; };

// shorthands
static inline Real sq(const Real &arg);
//{ return arg * arg; }

static inline Real cosh(const Real &arg)
{ Real e(exp(arg)); return (e + recip(e)) / 2; }

static inline Real sinh(const Real &arg)
{ Real e(exp(arg)); return (e - recip(e)) / 2; }

static inline Real tanh(const Real &arg)
{ return Real(1.0) - Real(2.0) / (exp(arg * 2) + 1); }

// comparison operators, semi-decidable
static inline bool operator < (const Real &lhs, const Real &rhs)
{ return (lhs - rhs).IsNegative(); }

static inline bool operator > (const Real &lhs, const Real &rhs)
{ return (lhs - rhs).IsPositive(); }

static inline bool operator != (const Real &lhs, const Real &rhs)
{ return (lhs - rhs).IsNonZero(); }

static inline char* NonZeroRealAsDecimal(const Real &arg, char *buffer, unsigned lenwanted) // force non-zero, then convert
{ if (arg.IsNonZero()) return arg.AsDecimal(buffer, lenwanted); 
  else return "n/a";    }

} // namespace

#endif
