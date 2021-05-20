CC=g++ 
CFLAGS= -std=gnu++11 -O3 -Wall -openmp
CXXFLAGS= -std=gnu++11 -O3 -Wall -openmp
CFILES=RunEmpirical.cpp mkmodel.cpp h_random.cpp graphmodel.cpp mcsim.cpp ioroutines.cpp optimizer.cpp
OBJ=RunEmpirical.o mkmodel.o h_random.o graphmodel.o mcsim.o ioroutines.o optimizer.o
DEPS=defs.hpp graphmodel.hpp h_random.hpp ioroutines.hpp mcsim.hpp mkmodel.hpp optimizer.hpp stat_tester.hpp

# every .cc file into a .o file
%.o : %.cc $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

runempirical: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)

clean:
	rm *.o
