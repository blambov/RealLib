CXX = g++
CXXFLAGS := -O3 

ifeq ($(USE_SSE),1)
  CXXFLAGS := $(CXXFLAGS) -msse2
endif

FLAGS_TO_PASS = \
	"CXX = $(CXX)" \
	"CXXFLAGS = $(CXXFLAGS)"

OBJECTS = convolution.o kernels.o DataManager.o LongFloat.o ErrorEstimate.o \
 RealEstimate.o RealFuncs.o RealObject.o Real.o MachineEstimate.o

all: Real.a makeexamples makemanual

Real.a: $(OBJECTS)
	ar -rvs Real.a $(OBJECTS)

MachineEstimate.h: MachineEstimateMul.h MachineEstimateSSE2.h
	if test "$(USE_SSE)" = "1"; then MACHEST=MachineEstimateSSE2; \
		else MACHEST=MachineEstimateMul; \
	fi; \
	cp $$MACHEST.h MachineEstimate.h

MachineEstimate.cpp: MachineEstimate.h MachineEstimateMul.cpp \
  MachineEstimateSSE2.cpp
	if test "$(USE_SSE)" = "1"; then MACHEST=MachineEstimateSSE2; \
		else MACHEST=MachineEstimateMul; \
	fi; \
	cp $$MACHEST.cpp MachineEstimate.cpp
	
makeexamples:
	cd examples; $(MAKE) $(FLAGS_TO_PASS)
	
makemanual:
	cd manual; $(MAKE) $(FLAGS_TO_PASS)
	
clean: 
	cd examples; $(MAKE) clean
	cd manual; $(MAKE) clean
	rm -f *.o *.a
	rm MachineEstimate.cpp
	rm MachineEstimate.h

convolution.o: convolution.cpp defs.h GCChelper.h convolution.h

kernels.o: kernels.cpp kernels.h defs.h GCChelper.h convolution.h

DataManager.o: DataManager.cpp DataManager.h defs.h GCChelper.h

ErrorEstimate.o: ErrorEstimate.cpp defs.h GCChelper.h ErrorEstimate.h \
  LongFloat.h

LongFloat.o: LongFloat.cpp LongFloat.h defs.h GCChelper.h DataManager.h \
  kernels.h

RealEstimate.o: RealEstimate.cpp defs.h GCChelper.h RealEstimate.h \
  LongFloat.h ErrorEstimate.h DataManager.h

RealFuncs.o: RealFuncs.cpp MachineEstimate.h defs.h GCChelper.h \
  RealEstimate.h LongFloat.h ErrorEstimate.h RealFuncs.h

MachineEstimate.o: MachineEstimate.cpp MachineEstimate.h defs.h \
  GCChelper.h RealEstimate.h LongFloat.h ErrorEstimate.h

RealObject.o: RealObject.cpp defs.h GCChelper.h RealObject.h \
  RealEncapsulation.h LongFloat.h RealEstimate.h ErrorEstimate.h \
  MachineEstimate.h RealFuncs.h

Real.o: Real.cpp defs.h GCChelper.h Real.h RealEncapsulation.h \
  LongFloat.h RealEstimate.h ErrorEstimate.h MachineEstimate.h \
  RealFuncs.h RealObject.h

