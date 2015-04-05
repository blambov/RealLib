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
#include "convolution.h"
#include <math.h>
#include <assert.h>

//using namespace std;

#define realindex(i)    ((i)*2)
#define imagindex(i)    ((i)*2+1)

// forward real-to-complex in-place fft
template <class TYPE>
void fft_fwd_ip_rc(int size, TYPE *weights, int *br, TYPE *data, int wstride);

// inverse complex-to-real in-place fft
template <class TYPE>
void fft_inv_ip_cr(int size, TYPE *weights, int *br, TYPE *data, int wstride);

// complex multiplication of vector
template <class TYPE>
void mul_complex(int size, TYPE *a, TYPE *b, double scale);


// returns the reversed bits of z for size = 2 ^ bits
static int bitreverse(int z, int bits, int size)
{
    int r = 0;
    for (int i=0; i<bits; ++i) {
        size >>= 1;
        if (z & 1) r += size;
        z >>= 1;
    }
    return r;
}

// constructor
template <class TYPE>
Convolution<TYPE>::Convolution(int size, TYPE PI2)
: m_Size(size)
{
    // verify size is correct. should be a power of 2
    assert ((size & (size - 1)) == 0);  
    int bits = int(log(double(size)) / log(2.0) + 0.5);

    // allocate memory
    m_pWeights = new TYPE [size*2];
    m_pBR = new int [size/2];
    assert(m_pWeights && m_pBR);

    // fill in weights and bit reverse vector
    int i;
    for (i=0; i<size/2; ++i) {
        m_pWeights[realindex(i)] = cos(PI2 * i / size);
        m_pWeights[imagindex(i)] = -sin(PI2 * i / size);

        m_pBR[i] = bitreverse(i, bits-1, size/2);
    }
}

// destructor
template <class TYPE>
Convolution<TYPE>::~Convolution()
{
    delete [] m_pWeights;
    delete [] m_pBR;
}

#ifndef __RESTRICT
#define restrict 
#endif

// convolve performs the operation.
// output goes to a. b is destroyed
template <class TYPE>
void Convolution<TYPE>::Convolve(TYPE * restrict a, TYPE * restrict b, int size)
{
    if (size==0) size = m_Size;
    // forward ffts. remember rc multiplies both by additional factor of 2
    int wstride = m_Size / size;

    fft_fwd_ip_rc(size, m_pWeights, m_pBR, a, wstride);
    fft_fwd_ip_rc(size, m_pWeights, m_pBR, b, wstride);

    // DC and Nyguest share one complex value
    // should be multiplied separately
    a[0] *= b[0] / (size * 4);
    a[1] *= b[1] / (size * 4);

    mul_complex(size/2-1, a+2, b+2, 1.0 / (size * 4));

    // inverse fft
    fft_inv_ip_cr(size, m_pWeights, m_pBR, a, wstride);
}

// multiply two complex vectors, in-place
template <class TYPE>
void mul_complex(int size, TYPE * restrict a, TYPE * restrict b, double scale)
{
    for (int i=0; i<size; ++i) {
        TYPE re(a[realindex(i)]*b[realindex(i)] - a[imagindex(i)]*b[imagindex(i)]);
        TYPE im(a[realindex(i)]*b[imagindex(i)] + a[imagindex(i)]*b[realindex(i)]);
        a[realindex(i)] = scale * re;
        a[imagindex(i)] = scale * im;
    }
}


// gentleman-sande decimation-in-frequency forward in-place fft
template <class TYPE>
void fft_fwd_ip(int size, TYPE * restrict weights, TYPE * restrict a, int wstride)
{
    for (int L=size; L>1; L >>= 1) {
        int r = size / L;
        int L2 = L >> 1;

        for (int j=0; j<L2; ++j) {
            TYPE wr = weights[realindex(j * r * wstride)];
            TYPE wi = weights[imagindex(j * r * wstride)];

            for (int k=0; k<r; ++k) {
                TYPE cr = a[realindex(k * L + j)];
                TYPE ci = a[imagindex(k * L + j)];
                TYPE dr = a[realindex(k * L + L2 + j)];
                TYPE di = a[imagindex(k * L + L2 + j)];

                a[realindex(k * L + j)] = cr + dr;
                a[imagindex(k * L + j)] = ci + di;

                cr -= dr;
                ci -= di;

                a[realindex(k * L + L2 + j)] = wr * cr - wi * ci;
                a[imagindex(k * L + L2 + j)] = wr * ci + wi * cr;
            }
        }
    }

    // permutation is skipped
}

// cooley-tukey decimation-in-time inverse in-place fft
template <class TYPE>
void fft_inv_ip(int size, TYPE * restrict weights, TYPE *restrict a, int wstride)
{
    // permutation is skipped

    for (int L = 2; L <= size; L <<= 1) {
        int r = size / L;
        int L2 = L / 2;

        for (int j=0; j<L2; ++j) {
            TYPE wr = weights[realindex(j * r * wstride)];
            TYPE wi = -weights[imagindex(j * r * wstride)]; // inverse

            for (int k=0; k<r; ++k) {
                TYPE cr = a[realindex(k * L + j)];
                TYPE ci = a[imagindex(k * L + j)];
                TYPE dr = a[realindex(k * L + L2 + j)];
                TYPE di = a[imagindex(k * L + L2 + j)];
                TYPE tr = wr * dr - wi * di;
                TYPE ti = wr * di + wi * dr;

                a[realindex(k * L + j)] = cr + tr;
                a[imagindex(k * L + j)] = ci + ti;

                a[realindex(k * L + L2 + j)] = cr - tr;
                a[imagindex(k * L + L2 + j)] = ci - ti;
            }
        }
    }
}

// real-to-complex step after fft_fwd
// the result is multiplied by 2
template <class TYPE>
void fft_realtocomplex(int size, TYPE * restrict weights, int * restrict br, TYPE * restrict a, int wstride)
{
    int size2 = size/2;
    TYPE pr, pi, mr, mi;
    int i, j;

    // calculate DC and Nyguest (the value at the center frequency)
    // both are real numbers, to avoid needing extra space they share one
    // complex point
    pr = a[realindex(0)]; pi = a[imagindex(0)];
    a[realindex(0)] = (pr + pi) * 2.0;
    a[imagindex(0)] = (pr - pi) * 2.0;

    // this is in the middle, bitreverse(size/2) == 1
    mr = a[realindex(1)]; mi = a[imagindex(1)];
    a[realindex(1)] = mr * 2.0;
    a[imagindex(1)] = mi * -2.0;

    // from here on, indexes are retrieved bitreversed
    // br(i*wstride) is the proper br(i) when the size is divided by wstride
    for (i=wstride, j=(size-1)*wstride; i<size2*wstride; (i+=wstride), (j-=wstride)) {
        pr = a[realindex(br[i])] + a[realindex(br[j])];
        pi = a[imagindex(br[i])] + a[imagindex(br[j])];
        mr = a[realindex(br[i])] - a[realindex(br[j])];
        mi = a[imagindex(br[i])] - a[imagindex(br[j])];

        a[realindex(br[i])] = pr + weights[realindex(i)] * pi + weights[imagindex(i)] * mr;
        a[imagindex(br[i])] = mi - weights[realindex(i)] * mr + weights[imagindex(i)] * pi;
        a[realindex(br[j])] = pr - weights[realindex(i)] * pi - weights[imagindex(i)] * mr;
        a[imagindex(br[j])] = -mi - weights[realindex(i)] * mr + weights[imagindex(i)] * pi;
    }
}

// complex-to-real step before fft_inv
template <class TYPE>
void fft_complextoreal(int size, TYPE * restrict weights, int * restrict br, TYPE * restrict a, int wstride)
{
    int size2 = size/2;
    TYPE pr, pi, mr, mi, zr, zi;
    int i, j;

    // DC and Nyquest were calculated using a different formula
    pr = a[realindex(0)]; pi = a[imagindex(0)];
    a[realindex(0)] = (pr + pi);
    a[imagindex(0)] = (pr - pi);

    // this is in the middle, bitreverse(size/2) == 1
    mr = a[realindex(1)]; 
    mi = a[imagindex(1)];
    a[realindex(1)] = mr * 2.0;
    a[imagindex(1)] = mi * -2.0;

    // from here on, indexes are retrieved bitreversed
    for (i=wstride, j=(size-1)*wstride; i<size2*wstride; (i+=wstride), (j-=wstride)) {
        pr = a[realindex(br[i])] + a[realindex(br[j])];
        pi = a[imagindex(br[i])] - a[imagindex(br[j])];
        mi = a[realindex(br[i])] - a[realindex(br[j])];
        mr = a[imagindex(br[i])] + a[imagindex(br[j])];

        zr = mr * weights[realindex(i)] - mi * weights[imagindex(i)];
        zi = mi * weights[realindex(i)] + mr * weights[imagindex(i)];

        a[realindex(br[i])] = pr - zr;
        a[imagindex(br[i])] = pi + zi;
        a[realindex(br[j])] = pr + zr;
        a[imagindex(br[j])] = zi - pi;
    }
}

template <class TYPE>
void fft_fwd_ip_rc(int size, TYPE * restrict weights, int * restrict br, TYPE * restrict a, int wstride)
{
    // perform a complex-to-complex fft on the data
    fft_fwd_ip(size / 2, weights, a, 2*wstride);

    // then use an additional step to get the actual result
    fft_realtocomplex(size / 2, weights, br, a, wstride);
}

template <class TYPE>
void fft_inv_ip_cr(int size, TYPE * restrict weights, int * restrict br, TYPE * restrict a, int wstride)
{
    // revert the operation of fft_realtocomplex
    fft_complextoreal(size / 2, weights, br, a, wstride);

    // perform a complex-to-complex fft
    fft_inv_ip(size / 2, weights, a, 2*wstride);
}


// instantiate
template
class Convolution<double>;


/*
test code

#include <stdio.h>

void main()
{
#define size 128
#define br(x) conv.m_pBR[x]
    Convolution<double> conv(size*2);
    double a[size*2];
    int i;

    for (i=0; i<size; ++i) {
        a[realindex(i)] = i + 1;
        a[imagindex(i)] = -i;
    }

    for (i=0; i<size; ++i) {
        printf("%3x: %lf+i%lf\n", i, a[realindex(i)], a[imagindex(i)]);
    }
    getchar();
    for (i=0; i<size; ++i) {
        printf("weights %3x: %lf+i%lf\n", i, conv.m_pWeights[realindex(i)], conv.m_pWeights[imagindex(i)]);
    }
    getchar();

    fft_fwd_ip(size, conv.m_pWeights, a, 2);

    for (i=0; i<size; ++i) {
        printf("%3x, %3x: %lf+i%lf\n", i, br(i), a[realindex(br(i))], a[imagindex(br(i))]);
    }
    getchar();

    fft_inv_ip(size, conv.m_pWeights, a, 2);

    for (i=0; i<size; ++i) {
        printf("%3x: %lf+i%lf\n", i, a[realindex(i)]/size, a[imagindex(i)]/size);
    }
    getchar();

    int sizer = size*2;
    for (i=0; i<sizer; ++i) a[i] = i;
    for (i=0; i<sizer; ++i) printf("%3x: %lf\n", i, a[i]);  
    getchar();

    fft_fwd_ip_rc(sizer, conv.m_pWeights, conv.m_pBR, a);

    for (i=0; i<size; ++i) {
        printf("%3x, %3x: %lf+i%lf\n", i, br(i), a[realindex(br(i))], a[imagindex(br(i))]);
    }
    getchar();

    fft_inv_ip_cr(sizer, conv.m_pWeights, conv.m_pBR, a);

    for (i=0; i<sizer; ++i) printf("%3x: %lf\n", i, a[i]/(sizer*2));    
    getchar();

    double b[size*2];
    b[0] = 1; b[1] = 1;
    for (i=2;i<sizer;++i) b[i] = 0;

    conv.Convolve(a, b);
    for (i=0; i<sizer; ++i) printf("%3x: %lf\n", i, a[i]/(sizer*2));    

    getchar();
}

 */

