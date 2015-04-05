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

#ifndef FILE_REAL_EXCEPTIONS_H
#define FILE_REAL_EXCEPTIONS_H

namespace RealLib {

class RealLibException : public std::exception {
    char m_what[128];
public:
    RealLibException(const char *what = NULL) throw();
    virtual const char *what() const throw()
              {  return m_what; }
    virtual const char *kind() const throw()
            {  return "RealLibException"; }
};

class PrecisionException : public RealLibException {
public:
    PrecisionException(const char *what = NULL) throw()
    : RealLibException(what) {}
    const char *kind() const throw()
            {  return "PrecisionException"; }
};

class DomainException : public RealLibException {
public:
    DomainException(const char *what = NULL) throw()
    : RealLibException(what) {}
    virtual const char *kind() const throw()
            {  return "DomainException"; }
};

static inline
std::ostream & operator << (std::ostream &os, RealLibException &e)
{
    return os << e.kind() << ": " << e.what();
}

}

#endif // FILE
