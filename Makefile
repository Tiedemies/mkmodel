CC=g++ 
CFLAGS= -std=gnu++11 -O3 -Wall -pedantic -fopenmp -lboost_date_time
CXXFLAGS= -std=gnu++11 -O3 -Wall -pedantic -fopenmp -lboost_date_time
CFILES=RunEmpirical.cpp mkmodel.cpp h_random.cpp graphmodel.cpp mcsim.cpp ioroutines.cpp optimizer.cpp test_io.cpp defs.cpp
EMPOBJ=RunEmpirical.o mkmodel.o h_random.o graphmodel.o mcsim.o ioroutines.o optimizer.o defs.o
TIOBJ=defs.o test_io.o ioroutines.o stat_tester.o graphmodel.o
STAT_P=ioroutines.o defs.o stat_tester.o graphmodel.o
DEPS=defs.hpp graphmodel.hpp h_random.hpp ioroutines.hpp mcsim.hpp mkmodel.hpp optimizer.hpp stat_tester.hpp

# every .cc file into a .o file
%.o : %.cc $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

RunEmpirical: $(EMPOBJ)
	$(CC) -o $@ $^ $(CFLAGS)

test_io: $(TIOBJ)
	$(CC) -o $@ $^ $(CFLAGS)

clean:
	rm *.o
