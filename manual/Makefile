OBJECTS = fragments.o harm.o harmfun.o \
	exp.o expimpr.o empty.o

EXAMPLES = fragments.out harm.out harmfun.out \
	exp.out expimpr.out empty.out
	
all: $(OBJECTS) $(EXAMPLES)

noobjects: $(OBJECTS)

%.out: %.o ../Real.a
	$(CXX) $(CXXFLAGS) -o $@ $< ../Real.a -lm

%.o: %.cpp ../Real.h ../RealEstimate.h ../MachineEstimate.h
	$(CXX) $(CXXFLAGS) -c $< -I ..
	
clean:
	rm -f *.o *.out
