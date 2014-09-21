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

#include "kernels.h"
#include "convolution.h"

#include <assert.h>
#include <math.h>

namespace RealLib {

#ifndef max
static inline int max(int x, int y) { return x > y ? x : y; }
#endif

i32 g_KernelsPrecision = 0;

#ifdef MULTIPLY_BY_CONVOLUTION
#define PI2 (2*3.14159265358979323846264338327950288419)
// buffers for convolution
double *g_pConvBufA = 0;
double *g_pConvBufB = 0;
// and the object itself
Convolution<double> *g_pConv = 0;
#endif

void InitializeKernels(u32 precision)
{
    g_KernelsPrecision = precision;

    #ifdef MULTIPLY_BY_CONVOLUTION
   if (precision >= CONVOLUTION_THRESHOLD) {
    u32 prec2pow;       
    // convolution will take 16 bits for a double
    // size should be at least twice the number of 16-bit words
    for (prec2pow = 16; prec2pow < precision*4; prec2pow <<= 1) ;
    
    g_pConvBufA = new double[prec2pow];
    g_pConvBufB = new double[prec2pow];
    g_pConv = new Convolution<double>(prec2pow, PI2);
   }
    #endif
}

void FinalizeKernels()
{
    #ifdef MULTIPLY_BY_CONVOLUTION
    if (g_KernelsPrecision >= CONVOLUTION_THRESHOLD) {
    delete g_pConvBufA;
    g_pConvBufA = 0;
    delete g_pConvBufB;
    g_pConvBufB = 0;
    delete g_pConv;
    g_pConv = 0;
    }
    #endif

    g_KernelsPrecision = 0;
}

// NormalizeMantissa: substraction can introduce zeroes in most
// significant positions of the mantissa. This function corrects
// such mantissas and returns the value that has to be substracted
// from the exponent. If this value is equal to working precision,
// the substraction function must recognize the value as Zero.
u32 NormalizeMantissa(u32 *man)
{
    i32 prec = g_KernelsPrecision - 1;
    if (man[prec] != 0) return 0;

    // find the first non-zero word
    while (--prec >= 0 && man[prec] == 0) ;
    
    // calculate needed offset
    prec = g_KernelsPrecision - (prec + 1);

    // do we have something to save?
    if (prec != g_KernelsPrecision) {
        i32 u;

        for (u=g_KernelsPrecision - 1; u >= prec; --u)
            man[u] = man[u - prec];

        for (;u >= 0; --u)
            man[u] = 0;
    }

    return prec;
}

// perform the actual addition
// returns true if there is carry
bool AddMantissa(u32 *man,          // destination, pre-initialized by a call to the default constructor
                 const u32 *full,   // the greater value
                 const u32 *part,   // the partial value,
                 i32 start)             // which is shifted by this many words
{
    int carry = 0;

    // start with carry if highest bit in what's left out is 1
    if (start != 0 && start <= g_KernelsPrecision)
        carry = part[start-1] >= (1u << 31);

    u64 v; i32 u;

    // add words
    for (u=0; u<g_KernelsPrecision - start; ++u) {
        v = u64(full[u]) + u64(part[u + start]) + u64(carry);
        man[u] = u32(v);
        carry = u32(v >> 32);
    }

    // update for carry
    for (; carry && u<g_KernelsPrecision; ++u) {
        man[u] = full[u] + carry;
        carry = man[u] == 0;
    }

    // just copy
    for (; u<g_KernelsPrecision; ++u) {
        man[u] = full[u];
    }

    return !!carry;
}

// adjust for calculations that don't fit the preallocated space.
// an extra pass might be needed if the leftover word introduces more carry.
// returns number of shifts done
u32  AdjustForCarry(u32 *man, 
                    u32 msw)            // most significant word, the one that doesn't fit in
{
    // round what's left over
    u32 carry = man[0] >= (1u<<31);
    i32 u;

    // shift
    for (u = 1; u<g_KernelsPrecision && carry; ++u) {
        man[u-1] = man[u] + 1;
        carry = man[u-1] == 0;
    }

    for (; u<g_KernelsPrecision; ++u)
        man[u-1] = man[u];

    // put new value
    man[u-1] = msw + carry;

    // reiterate if necessary
    if (man[u-1] == 0) return 1 + AdjustForCarry(man, 1);
    else return 1;
}

// perform the actual substraction
// returns true if part was greater and the result must be negated
bool SubMantissa(u32 *man,          // destination, pre-initialized by a call to the default constructor
                 const u32 *full,   // the greater value
                 const u32 *part,   // the partial value,
                 i32 start)             // which is shifted by this many words
{
    int carry = 0;

    // start with carry if highest bit in what's left out is 1
    if (start != 0 && start <= g_KernelsPrecision)
        carry = part[start-1] >= (1u << 31);

    u64 v; i32 u;

    // sub words
    for (u=0; u<g_KernelsPrecision - start; ++u) {
        v = u64(full[u]) - u64(part[u + start]) - u64(carry);
        man[u] = u32(v);
        carry = (v >> 32) != 0;
    }

    // update for carry
    for (; carry && u<g_KernelsPrecision; ++u) {
        man[u] = full[u] - carry;
        carry = man[u] == u32(-1);
    }

    // just copy
    for (; u<g_KernelsPrecision; ++u) {
        man[u] = full[u];
    }

    return !!carry;
}

// negate a mantissa. needed if SubMantissa returned true
void NegMantissa(u32 *man)
{
    for (i32 u = 0; u < g_KernelsPrecision; ++u)
        man[u] = ~man[u];
}

// MulMantissa
// perform actual multiplication
// the most significant word of the result is not put in man.
// instead it is returned, so no precision will be lost if
// it is zero
#ifndef MULTIPLY_BY_CONVOLUTION
u32 MulMantissa(u32 *man, 
                const u32 *a,
                const u32 *b,
                i32 inputstart,
                i32 inputlen)
#else
u32 MulMantissaDirect(u32 *man, 
                const u32 *a,
                const u32 *b,
                i32 inputstart,
                i32 inputlen)
#endif
{
    u32 carry = 0;
    u64 u, w = 0;
    i32 i = g_KernelsPrecision - inputlen * 2 + 1;
    i32 j;
    i32 k = 0;

    // start by only calculating carry
    for (; i < 0 && k < inputlen; ++i, ++k) {
        w >>= 32;
        w += u64(carry)<<32;
        carry = 0;
        for (j=0; j<=k; ++j) {
            u = u64(a[j+inputstart]) * u64(b[k-j+inputstart]);
            w += u;
            if (w < u) ++carry;     // this is a trick to check for carry in u64s
                                    // assembler would've made this a lot easier
        }
    }
    // alternatively
    for (j=0; j<i; ++j) man[j] = 0;
    
    assert(i>=0);

    // we didn't write till now. 
    // besides carry, we should add 1 if the previous value had 1 in MS bit
    if (w & 0x80000000) w += u64(1)<<32;

    // start writing
    for (; k < inputlen; ++i, ++k) {
        w >>= 32;
        w += u64(carry)<<32;
        carry = 0;
        for (j=0; j<=k; ++j) {
            u = u64(a[j+inputstart]) * u64(b[k-j+inputstart]);
            w += u;
            if (w < u) ++carry;
        }
        man[i] = u32(w);
    }

    for (; i < g_KernelsPrecision; ++i, ++k) {
        w >>= 32;
        w += u64(carry)<<32;
        carry = 0;
        for (j=k - inputlen + 1; j<inputlen; ++j) {
            u = u64(a[j+inputstart]) * u64(b[k-j+inputstart]);
            w += u;
            if (w < u) ++carry;
        }
        man[i] = u32(w);
    }
    w >>= 32;
    assert(!carry);

    // leave the last word as return value
    return u32(w);
}

#ifdef MULTIPLY_BY_CONVOLUTION
u32 MulMantissa(u32 *man, 
                const u32 *a,
                const u32 *b,
                i32 inputstart,
                i32 inputlen)
{
    // do it directly if it would be faster
    if (inputlen < CONVOLUTION_THRESHOLD)       
        return MulMantissaDirect(man, a, b, inputstart, inputlen);

    int i;
    int prec = inputlen;
    int prec2pow = 16;
    if (inputlen == g_KernelsPrecision)
       prec2pow = g_pConv->GetSize();
    else while (prec2pow < prec*4) prec2pow *= 2;

    // initialize buffers to input
    for (i=0; i<inputlen; ++i) {
        g_pConvBufA[i*2] = a[i+inputstart] & 0xFFFF;
        g_pConvBufA[i*2+1] = (a[i+inputstart] >> 16) & 0xFFFF;
    }
    i*=2;
    while (i<prec2pow) g_pConvBufA[i++] = 0;

    for (i=0; i<inputlen; ++i) {
        g_pConvBufB[i*2] = b[i+inputstart] & 0xFFFF;
        g_pConvBufB[i*2+1] = (b[i+inputstart] >> 16) & 0xFFFF;
    }
    i*=2;
    while (i<prec2pow) g_pConvBufB[i++] = 0;

    // convolve
    g_pConv->Convolve(g_pConvBufA, g_pConvBufB, prec2pow);

    // make each value 16-bit
    double carry = 0, t;
    for (i=0; i<inputlen - 1-inputstart; ++i) {
        t = floor(g_pConvBufA[i] + carry + 0.5);    // round it too
        carry = floor(ldexp(t, -16));
        g_pConvBufA[i] = t - ldexp(carry, 16);
    }
    // from here on we start writing, one in MSB of previous word is carry
    if (g_pConvBufA[i-1] > (1<<15)) carry += 1;

    for (; i<(prec + inputlen) * 2; ++i) {
        t = floor(g_pConvBufA[i] + carry + 0.5);    // round it too
        carry = floor(ldexp(t, -16));
        g_pConvBufA[i] = t - ldexp(carry, 16);
    }

   for (i=0;i<=inputstart-inputlen;++i) man[i] = 0;
    // write the result
    for (i=max(0, inputlen-1-inputstart); i<prec + inputlen - 1; ++i) {
        man[i-inputlen+1+inputstart] = u32(g_pConvBufA[i*2]) + (u32(g_pConvBufA[i*2+1]) << 16);
    }

    // leave the last word out
    return u32(g_pConvBufA[i*2]) + (u32(g_pConvBufA[i*2+1]) << 16);
}
#endif

bool MultipliedByConvolution(i32 inputlen)
{
   return inputlen >= CONVOLUTION_THRESHOLD;
}

template<class T>
static inline void swap(T &a, T &b) { T c=a;a=b;b=c; }

// auxilliary function to help division
// amsw is the most significant word of a
// aofs is how many words a is shifted, with the msw's
// taken as 0, the first substituted by the amsw
// the result is shifted aofs words
// bscale is assumed positive, < 32
// returns false is a was < b, possibly breaking with an
// incomplete res.
bool SubManBScaled(u32 *res,
                   const u32 *a,
                   const u32 *b,
                   u32 &amsw,
                   i32 bscale,
                   i32 inputlen,
                   i32 inputstart,
                   i32 aofs)
{
#define combinewords(a, b, bscale) (bscale==0 ? b : ((a >> (32 - bscale)) + (b << bscale)))
    int carry = 0;

    u64 v; 
    i32 u; 
    u32 s = combinewords(0, (b[inputstart]), bscale);
    
    for (u=inputstart; u<inputstart+aofs; ++u) {
        v = u64(0) - u64(s) - u64(carry);
        res[u] = u32(v);
        carry = (v >> 32) != 0;
        s = combinewords(b[u], b[u+1], bscale);
    }       
    
    
    // sub words
    for (u=0; u<inputlen-1; ++u) {
        v = u64(a[u-aofs+inputstart]) - u64(s) - u64(carry);
        res[u+inputstart] = u32(v);
        carry = (v >> 32) != 0;
        s = combinewords(b[inputstart+u], b[inputstart+u+1], bscale);
    }
    
   {
        v = u64(a[u-aofs+inputstart]) - u64(s) - u64(carry);
        res[u+inputstart] = u32(v);
        carry = (v >> 32) != 0;
        s = combinewords(b[inputstart+u], 0, bscale);
    }
    
    v = u64(amsw) - u64(s) - u64(carry);
    carry = (v >> 32) != 0;
    
    if (carry) return false;
    else {
       amsw = u32(v);
       return true;
    }
}

i32 DivMantissa(u32 *man,
                const u32 *a,
                const u32 *b,
                i32 inputstart,
                i32 inputlen,
                u32 *temp1,
                u32 *temp2)
{
   u32 amsw = 0;
   i32 sc;
   u32 r = 0;
   i32 e = 1;
   i32 j = inputstart + inputlen - 1;
   i32 i = j;
   i32 ofs = 0;
   
   for (int k=0;k<inputstart;++k) man[k] = 0;
   
   for (sc=31; sc>=0; --sc)
      if (SubManBScaled(temp1, a, b, amsw, sc, inputlen, inputstart, 0)) break;
   
   if (sc<0) {
      e = 0; --i;
      amsw = a[inputlen-1 + inputstart];
      for (sc = 31; sc>=0; --sc)
         if (SubManBScaled(temp1, a, b, amsw, sc, inputlen, inputstart, 1)) break;      

      assert(sc >= 0);
   }
   
   r |= 1<<sc;
   while (j>=inputstart) {
      while (--sc >= 0) {
         if (SubManBScaled(temp2, temp1, b, amsw, sc, inputlen, inputstart, 0)) {
            r |= 1<<sc;
            swap(temp1, temp2);
         }
      }
      
      ofs = 0;
      while (sc < 0 && j >= inputstart) {
         ++ofs;
         sc = 32; --i;
         man[j--] = r;
         if (j<inputstart) break;
         
         amsw = temp1[inputlen - ofs + inputstart];
         r = 0;
         
         for (sc = 31; sc>=0; --sc)
               if (SubManBScaled(temp2, temp1, b, amsw, sc, inputlen, inputstart, ofs)) {
               r |= 1<<sc;
               swap(temp1, temp2);
               break;      
            }
      }
   }
   
   // check if we need to round up
   if (SubManBScaled(temp2, temp1, b, amsw, 31, inputlen, inputstart, ofs))
   {
      while (++j < g_KernelsPrecision && ++man[j] == 0) ;
      
      if (j==g_KernelsPrecision) {  // carry on msw means we have 1(0)
         ++e;
         man[j-1] = 1;
      }
   }
   
   return e;
}

// scale mantissa: multiplication by u32 multiplier
// implemented for performance
u32 ScaleMantissa(u32 *man,
                  const u32 *src,
                  u32 multiplier)
{
    u64 v = 0;

    for (i32 i=0; i<g_KernelsPrecision; ++i)
    {
        v += u64(src[i]) * u64(multiplier);
        man[i] = u32(v);
        v >>= 32;
    }

    return u32(v);
}

i32 InvScaleMantissa(u32 *man,
                     const u32 *src,
                     u32 divisor)
{
   i32 i = g_KernelsPrecision - 1;
   i32 j = i;
   i32 e = 0;
   
   u64 v = src[i];
   if (v < divisor) {
      v = (v << 32) + src[--i];
      e = -1;
   }
   
   while (i > 0)
   {
      man[j--] = u32(v / divisor);
      v = ((v % divisor) << 32) + src[--i];
   }
   
   man[j--] = u32(v / divisor);
   
   if (j == 0) {     // this would happen if msw in src was < divisor
      v = (v % divisor) << 32;
      man[j--] = u32(v / divisor);
   }
   
   // round the result; j is -1
   if ((v % divisor) > divisor/2) {
      while (++j < g_KernelsPrecision && ++man[j] == 0) ;
      
      if (j==g_KernelsPrecision) {  // carry on msw means we have 1(0)
         ++e;
         man[j-1] = 1;
      }
   }
   return e;
}

// binary scale mantissa, i.e. multiply by 1<<scale, where scale < 32
u32 BScaleMantissa(u32 *man,
                  const u32 *src,
                  u32 scale)
{
    u32 v = 0;

    for (i32 i=0; i<g_KernelsPrecision; ++i)
    {
        man[i] = (src[i] << scale) | v;
        v = src[i] >> (32 - scale);
    }

    return u32(v);
}

}
