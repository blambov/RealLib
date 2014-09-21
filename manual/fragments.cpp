#include <iostream>
#include <iomanip>
#include "Real.h"
using namespace std;
using namespace RealLib;

int main() {
	InitializeRealLib();
	Real a, b(4);
	b = Pi / b;
	a = sin(b);
	cout << "sin(Pi/4) is " << a << endl;
	// ...
	Real c(a * 2 / sqrt(Real(2))); 
	for (int i=0; i<1000; ++i)
		c = sin(c);
	cout << "sin(sin(...(1)...)) (1000 times) is " << c << endl;
	// ...
	cout << " or " << scientific << setprecision(120) << c << endl;
	// ...
	cout << "the double representation of 0.1 is " 
		<< fixed << noshowpoint << Real(0.1) << endl;
	cout << "and its distance from 0.1 is "  
		<< showpos << scientific << uppercase 
		<< showpoint << Real(0.1) - Real("0.1") << endl;
	cout << "double(0.1) < Real(\"0.1\"): " 
		<< boolalpha << (Real(0.1) < Real("0.1")) << endl;
	cout << "0.1 converted to double is "
		<< noshowpos << fixed << Real("0.1").AsDouble() << endl;
	// ...
	FinalizeRealLib();
	return 0;
}
