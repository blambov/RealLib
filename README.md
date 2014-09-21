RealLib: Library for exact real number computations
=======

RealLib is a real number computation package in C++. Its primary aim is to be efficient and to avoid the huge overheads usually associated with real number computations.

RealLib3, which is the current version of the library, makes this possible. In this version of the library, the user can work on both a layer where numbers are represented as terms describing the computation, and a layer where functions on real numbers compute on the level of approximations to the real number. The latter can be very fast, depending on the precision that is actually required by the computation. In the cases where machine precision is sufficient, exact real computations can be executed at a speed comparable to the speed of double precision arithmetic, sometimes even running on par with it.

Version 3 comes with a detailed user manual that describes the system from the user's point of view and shows how such fast computations are made possible. An introductory paper, which contains performance comparison to two other real number packages and to machine precision floating point, is also available. The theoretical background, an approach to use Type-1 functions for real number computations, is described in a separate paper.
