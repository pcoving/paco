CXX = g++44
CXXFLAGS = -O3 -fopenmp

paco: paco.o 
	$(CXX) $(CXXFLAGS) -o paco paco.o 	

paco.o: paco.cpp paco.hpp quad.hpp
	$(CXX) $(CXXFLAGS) -c paco.cpp	

clean:
	rm *.o