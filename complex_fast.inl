#ifndef __COMPLEX_FAST_INL
#define __COMPLEX_FAST_INL

namespace RealLib {

template <class T> class complex {
public:
      T re;
      T im;
      
      complex(const T re_=T(), const T im_=T()) : re(re_), im(im_) {}
      //complex(const complex<T> &s) : re(s.re), im(s.im) {}
      
	  template <class X> 
      explicit complex(const complex<X> &s) : re(s.re), im(s.im) {}
      
      T real() const {return re;} 
      T imag() const {return im;}
      
      inline complex<T> operator= (const T re_) {re=re_;im=T(0.0f);return *this;}
      inline complex<T> operator+=(const T re_) {re+=re_;return *this;}
      inline complex<T> operator-=(const T re_) {re-=re_;return *this;}
      inline complex<T> operator*=(const T re_) {re*=re_;im*=re_;return *this;}
      inline complex<T> operator/=(const T re_) {re/=re_;im/=re_;return *this;}
      
      template <class X> inline 
      complex<T> operator= (const complex<X> s) 
            {re = s.re; im = s.im; return *this;}
      
      template <class X> inline
      complex<T> operator+=(const complex<X> s) 
            {re += s.re; im += s.im;return *this;}
      
      template <class X> inline
      complex<T> operator-=(const complex<X> s) 
            {re -= s.re; im -= s.im; return *this;}
      
      template <class X> inline
      complex<T> operator*=(const complex<X> s) 
            {     T v = re*s.im;
                  T w = im*s.im;
                  re *= s.re;
                  im *= s.re;
                  re -= w;
                  im += v;
                  return *this;
                  }
      
      template <class X> inline
      complex<T> operator/= (const complex<X> s) 
            {     T v = re*s.im;
                  T w = im*s.im;
                  T u = s.re*s.re + s.im*s.im;
                  re *= s.re;
                  im *= s.re;
                  re = (re + w)/u;
                  im = (im - v)/u;
                  return *this;
                  }
      
      inline void mulj() { T t=-im;im=re;re=t; }
};

template <class T>
inline complex<T> operator+(const complex<T> a, const complex<T> b)
{
      return complex<T> (a.re + b.re, a.im + b.im);
}

template <class T>
inline complex<T> operator-(const complex<T> a, const complex<T> b)
{
      return complex<T> (a.re - b.re, a.im - b.im);
}

template <class T>
inline complex<T> operator-(const complex<T> a)
{
      return complex<T> (-a.re, -a.im);
}

template <class T>
inline complex<T> operator<<(const complex<T> a, const int b)
{
      return complex<T> (a.re << b, a.im << b);
}

template <class T>
inline complex<T> operator*(const complex<T> a, const complex<T> b)
{
      return complex<T> (a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re);
}

template <class T>
inline complex<T> operator*(const complex<T> a, const T b)
{
      return complex<T> (a.re * b, a.im * b);
}

template <class T> inline
complex<T> operator/ (const complex<T> a, const complex<T> b) 
{
	T v = a.re*b.im;
	T w = a.im*b.im;
	T u = b.re*b.re + b.im*b.im;
	T re = a.re * b.re;
	T im = a.im * b.re;
	return complex<T> ((re + w)/u, (im - v)/u);
}
      
template <class T>
inline complex<T> operator/(const complex<T> a, const T b)
{
      return complex<T> (a.re / b, a.im / b);
}

template <class T> 
inline complex<T> exp(const complex<T> s) 
{
      T exp_x = exp(s.re);
      
      return complex<T> (exp_x * cos(s.im), exp_x * sin(s.im));
}

template <class T>
inline T abs(const complex<T> s)
{
      return sqrt(s.re*s.re + s.im*s.im);
}

template <class T>
inline complex<T> mulj(const complex <T> s)
{
      return complex<T> (-s.im, s.re);
}

template <class T>
inline complex<T> submulj(const complex <T> a, const complex<T> b)
{
      return complex<T> (b.im - a.im, a.re - b.re);
}

template <class T>
inline complex<T> addconj(const complex <T> a, const complex<T> b)
{
      return complex<T> (a.re + b.re, a.im - b.im);
}

template <class T>
inline complex<T> conj(const complex <T> a)
{
      return complex<T> (a.re, -a.im);
}

template <class T>
inline complex<T> subconj(const complex <T> a, const complex<T> b)
{
      return complex<T> (a.re - b.re, a.im + b.im);
}

}

#endif // __COMPLEX_FAST_INL

