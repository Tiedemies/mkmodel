CC=g++ 
CFLAGS= -std=gnu++11 -O3 -Wall -fopenmp -lboost_date_time
CXXFLAGS= -std=gnu++11 -O3 -Wall -fopenmp -lboost_date_time
CFILES=RunEmpirical.cpp mkmodel.cpp h_random.cpp graphmodel.cpp mcsim.cpp ioroutines.cpp optimizer.cpp test_io.cpp
OBJ=RunEmpirical.o mkmodel.o h_random.o graphmodel.o mcsim.o ioroutines.o optimizer.o test_io.o
DEPS=defs.hpp graphmodel.hpp h_random.hpp ioroutines.hpp mcsim.hpp mkmodel.hpp optimizer.hpp stat_tester.hpp

# every .cc file into a .o file
%.o : %.cc $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

runempirical: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)

test_io:
	g++ -o test_io test_io.cpp ioroutines.cpp -Wall -fopenmp -std=gnu++11

clean:
	rm *.o
