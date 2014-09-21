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

  convolution.h

  This file defines the convolution operation used for multiplication.

*/

#ifndef FILE_CONVOLUTION_H
#define FILE_CONVOLUTION_H


// real convolution object
// size is the convolution size, that is 
// for kernel length N and signal length M, N + M - 1
// the inputs to Convolve are of length size, padded with 0.
// the output of Convolve is in a; b is destroyed.

// this object works only on size = 2 ^ k for an integer k


template <class TYPE>
class Convolution {
private:
    TYPE *m_pWeights;       // pre-computed weights vector (exp(jpi2 * i/size))
    int *m_pBR;             // pre-computed bit-reversed vector
    int m_Size;             // size of the operation

public:
    Convolution(int size, TYPE PI2);    
    ~Convolution();         

    // perform convolution on the given arrays
    // output in a, b destroyed
    void Convolve(TYPE *a, TYPE *b, int size = 0);
    
    int GetSize() 
    { return m_Size; }
};

#endif
