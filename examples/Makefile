OBJECTS = taylor.o tayest.o taydbl.o pidigits.o \
	complex.o fibonacci.o convolve.o benchops.o \
	iterated.o linear.o layers.o

EXAMPLES = taylor.out tayest.out taydbl.out pidigits.out \
	complex.out fibonacci.out convolve.out benchops.out \
	iterated.out linear.out layers.out
	
all: $(OBJECTS) $(EXAMPLES)

noobjects: $(OBJECTS)

%.out: %.o ../Real.a
	$(CXX) $(CXXFLAGS) -o $@ $< ../Real.a -lm

%.o: %.cpp ../Real.h ../RealEstimate.h ../MachineEstimate.h
	$(CXX) $(CXXFLAGS) -c $< -I ..
	
clean:
	rm -f *.o *.out
