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

#include <iomanip>
#include "MachineEstimate.h"

namespace RealLib {

double MachineEstimate::plusinf = 0.0;
double MachineEstimate::minusinf = 0.0;
double MachineEstimate::OnePEps = 1.0;
double MachineEstimate::OneMEps = 1.0;
double MachineEstimate::pi_low = 3.14159265358979323846264338328;
double MachineEstimate::pi_high = 0.0;
double MachineEstimate::ln2_low = 0.0;
double MachineEstimate::ln2_high = 0.0;

#ifdef _MSC_VER
#define nextafter _nextafter
#endif

void MachineEstimate::BeginComputation()
{
    double zero = 0.0;
    plusinf = 1.0/zero;
    minusinf = -plusinf;
    OnePEps = nextafter(1.0, plusinf);
    OneMEps = nextafter(1.0, minusinf);
    pi_high = nextafter(pi_low, plusinf);
    ln2_low = nextafter(std::log(2.0), minusinf);
    ln2_high = nextafter(ln2_low, plusinf);
}

void MachineEstimate::FinishComputation()
{
}

template <>
MachineEstimate ln2(unsigned int prec)
{
    return MachineEstimate(MachineEstimate::ln2_low, MachineEstimate::ln2_high);
}

template <>
MachineEstimate pi(unsigned int prec)
{
    return MachineEstimate(MachineEstimate::pi_low, MachineEstimate::pi_high);
}

std::ostream& operator <<(std::ostream &os, const MachineEstimate &e)
{   
    return os << e.weak_AsDouble();
    /*
    if (e.low < 0 && e.high > 0) return os << 0.0;
    else {
        double d(e.weak_AsDouble()); 
        if (d < 1) os << d;
        else {
            int op = os.precision();
            double v = d;
            int np = op;
            while (v >= 1 && np > 0) {
                --np;
                v /= 10;
            }
            os << std::setprecision(np) << d << std::setprecision(op);
        }
    }*/
}

}
