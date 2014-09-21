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
