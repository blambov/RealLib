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
