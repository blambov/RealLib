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

  RealObject.h

  This header defines the objects that represent real numbers.
  A real is a tree that represents the procedure by which it
  is calculated. It might be a conversion from another type
  (double or string), a constant (RealNullary), or a operation
  applied to another real (RealUnary and RealBinary).

  RealObjects save this information, also handle the repeated
  usage of a term, and the caching of calculated estimates.
  They are written in zero-overhead manner, i.e. care is taken
  to never introduce extra calculation or extra storage than
  what is minimally needed, and what is actually needed.

*/

#ifndef FILE_REALOBJECT_H
#define FILE_REALOBJECT_H

#include "RealEncapsulation.h"

namespace RealLib {

struct ObjectList;
struct EvalList;
    
// g_pEstimatesList stores pointers to all temporary estimates
// it is used to delete them in case of library Reset (precision
// change), when they become invalid
extern ObjectList *g_pEstimatesList;

typedef const char* (*OracleFunction) (unsigned precision);

// prec and the Encapsulation working precision will grow together
    
// RealObject: the abstract base class

class RealObject {
    // reference counter
    i32 m_RefCount;
    // a pointer to the pointer to this object in g_pEstimatesList
    ObjectList *m_pPtrInObjList;

protected:
    // temporary Encapsulation
    Encapsulation *m_pEstimate;

    RealObject(u32 rc = 0);
    
    // evaluation procedure. 
    // To be overridden by implementations
    virtual Encapsulation Evaluate() = 0;

    // destructor
    virtual ~RealObject();

    // release the other objects this one points to
    // only in RealUnary and RealBinary
    virtual void ReleaseSiblings();

public:
    // depth. equal to the max depth of the siblings + 1
    i32 m_Depth;
    // Encapsulation refs. set to refcount when the Encapsulation is created.
    // decreased by GetEstimate. The Encapsulation is destroyed when
    // m_EstimateRefs reaches 0.
    i32 m_EstimateRefs;

    // get the Encapsulation, evaluate if necessary
    Encapsulation GetEstimate();
   // do we hold an Encapsulation?
   bool HasEstimate()
   { return !!m_pEstimate; }

   void AddToEstimatesList();
   void RemoveFromEstimatesList();

    // create a new Encapsulation (differs only in arrays, usually new Encapsulation (val))
   // and add to g_pEstimatesList
   virtual Encapsulation *CreateEstimate(const Encapsulation &val);
   // destroy the Encapsulation and remove this object from g_pEstimatesList
    virtual void DestroyEstimate();

    // add reference. Called when a pointer to this object is saved
    RealObject* AddRef();
    // release: decreases refcount and deletes object if rc == 0
    void Release();
    // get the reference count
    u32 GetRefCount() 
    { return m_RefCount; }

    // returns the sibling with the given index
    // returns NULL if there are no more
    // to be used to do recursion with a separate stack
    virtual RealObject* GetSibling(int index);
    // non-recursive release, to be combined with GetSibling
    void NonRecursiveRelease();
    
    /* obsolete older version
    
    // gets all objects that need to be pre-evaluated
    // for the evaluation procedure not to dig deeper
    // than the required level. calls testers to see
    // if search can be stopped. won't get any list if 
    // pList is Null. returning from recursion calls
    // OnExit
    typedef bool (RealObject::*TesterFunc)(i32 Depth);
    typedef void (RealObject::*OnExitFunc)();

    virtual void GetDepthList(i32 Depth, TesterFunc Test, OnExitFunc OnExit, EvalList **pList, bool DiveIntoDeeperFirst);
    
   // tester for depth-restricted evaluation
   // no OnExit needed
   bool HasEstimateReference(i32 Depth);  // stops if m_EstimateRefs > 0, marks if already visited
                                          // followed by evaluation which cleans marks
   
   // tester/OnExit pair for depth-restricted release
   bool StartRelease(i32 Depth); // decreases m_RefCount, stops if > 0, or if Depth < m_Depth with ReleaseSiblings()
   void FinishRelease();         // cleans-up object without releasing siblings if m_RefCount == 0
    */

   /*
   // testers if depth search can be stopped
    // decreases refcount, returns true if > 0
    bool WontBeDeleted();
    // not used
    // stops if m_EstimateRefs > 0. 
    // increases m_EstimateRefs if not to mark visited objects
    bool HasEstimateReference();
    // on exit function for deletion
    // deletes the object if m_RefCount == 0
    void NonRecursiveRelease();
    // not used
    void SetEstimateReference();
   */
};

struct ObjectList {
    ObjectList *next;
    ObjectList *prev;
    RealObject *obj;
};

// this struct saves RealObjects to be visited
struct EvalList {
    EvalList *next;
    bool used;
    RealObject *obj;
};

// real from double. Handles conversion.

class RealFromDouble : public RealObject {
    double m_Value;
protected:
    virtual ~RealFromDouble();
    virtual Encapsulation Evaluate();

public:
    RealFromDouble(const double src);
};

// real from string. Handles conversion
class RealFromString : public RealObject {
    char *m_pString;
protected:
    virtual ~RealFromString();
    virtual Encapsulation Evaluate();
public:
    RealFromString(const char *src);
};

// real from oracle function
class RealFromOracle : public RealObject {
    OracleFunction m_pOracle;
protected:
    virtual ~RealFromOracle();
    virtual Encapsulation Evaluate();
public:
    RealFromOracle(OracleFunction oracle);
};

// RealNullary. constants
class RealNullary : public RealObject {
public:
    typedef Encapsulation (*FuncNullary) (unsigned int prec, UserInt otherData);
private:
    FuncNullary m_pFunc;
    UserInt m_iUserData;
protected:
    virtual ~RealNullary();
    virtual Encapsulation Evaluate();
public:
    RealNullary(FuncNullary pFunc, UserInt userData);
};

// unary functions
class RealUnary : public RealObject {
public:
    typedef Encapsulation (*FuncUnary) (const Encapsulation &arg, UserInt otherData);
private:
    FuncUnary m_pFunc;
    RealObject *m_pArg;
    UserInt m_iUserData;

protected:
    virtual ~RealUnary();
    virtual void ReleaseSiblings();
    virtual Encapsulation Evaluate();
public:
    //virtual void GetDepthList(i32 Depth, TesterFunc Test, OnExitFunc OnExit, EvalList **pList, bool DiveIntoDeeperFirst);
    RealUnary(FuncUnary pFunc, RealObject *pArg, UserInt userData);
    virtual RealObject* GetSibling(int index);
};

// binary functions
class RealBinary : public RealObject {
public:
    typedef Encapsulation (*FuncBinary) (const Encapsulation &left, const Encapsulation &right, UserInt otherData);
private:
    FuncBinary m_pFunc;
    RealObject *m_pLeft;
    RealObject *m_pRight;
    UserInt m_iUserData;

protected:
    virtual ~RealBinary();
    virtual void ReleaseSiblings();
    virtual Encapsulation Evaluate();
public:
    //virtual void GetDepthList(i32 Depth, TesterFunc Test, OnExitFunc OnExit, EvalList **pList, bool DiveIntoDeeperFirst);
    RealBinary(FuncBinary pFunc, RealObject *pLeft, RealObject *pRight, UserInt userData);
    virtual RealObject* GetSibling(int index);
};

/*
class RealBinaryOnInt : public RealObject {
public:
    typedef Encapsulation (*FuncBinaryOnInt) (const Encapsulation &arg, long iarg);
private:
    FuncBinaryOnInt m_pFunc;
    RealObject *m_pArg;
    long m_iArg;
    
protected:
    virtual ~RealBinaryOnInt();
    virtual void ReleaseSiblings();
    virtual Encapsulation Evaluate();
public:
    //virtual void GetDepthList(i32 Depth, TesterFunc Test, OnExitFunc OnExit, EvalList **pList, bool DiveIntoDeeperFirst);
    RealBinaryOnInt(FuncBinaryOnInt pFunc, RealObject *pArg, long iarg);        
    virtual RealObject* GetSibling(int index);
};*/

// functions on array
// this is a bit complicated.
// first a RealArray object is created that holds references to all
// realobject's in the passed array
// then, for each element in the array, a RealArrayElement is
// created that holds a reference to the array.
// Evaluation of the element is passed along as evaluation of the
// array, having set the request index to its own.
class RealArray : public RealObject {
public:
    typedef void (*FuncArray) (ArrayInterface<Encapsulation>&, UserInt otherData);
   //typedef void (*FuncArray) (Encapsulation* arr, unsigned int count, void* userdata);
private:
   FuncArray m_pFunc;
   RealObject **m_pArray;
   Encapsulation *m_pEstimateArray;
   u32 m_uCount;
   UserInt m_iUserData;
   i32 m_uRequestIndex;

protected:
   virtual ~RealArray();
   virtual void ReleaseSiblings();
   virtual Encapsulation Evaluate();
public:
   virtual Encapsulation *CreateEstimate(const Encapsulation &val);
   virtual void DestroyEstimate();
   void SetRequestIndex(i32 index);

   //virtual void GetDepthList(i32 Depth, TesterFunc Test, OnExitFunc OnExit, EvalList **pList, bool DiveIntoDeeperFirst);
   RealArray(FuncArray pFunc, RealObject **pArray, unsigned int count, UserInt user);   
    virtual RealObject* GetSibling(int index);
};

class RealArrayElement : public RealObject {
private:
   RealArray* m_pArray;
   u32 m_uMyIndex;
protected:
   virtual ~RealArrayElement();
   virtual void ReleaseSiblings();
   virtual Encapsulation Evaluate();
public:
   virtual Encapsulation *CreateEstimate(const Encapsulation &val);
   virtual void DestroyEstimate();
   //virtual void GetDepthList(i32 Depth, TesterFunc Test, OnExitFunc OnExit, EvalList **pList, bool DiveIntoDeeperFirst);
   RealArrayElement(RealArray *pArray, u32 myindex);    
    virtual RealObject* GetSibling(int index);
};

} // namespace

#endif
