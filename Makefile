runempirical: RunEmpirical.cpp mkmodel.cpp h_random.cpp graphmodel.cpp mcsim.cpp ioroutines.cpp optimizer.cpp
	g++ -std=gnu++11 -O3 -o RunEmpirical RunEmpirical.cpp mkmodel.cpp h_random.cpp graphmodel.cpp mcsim.cpp ioroutines.cpp optimizer.cpp -fopenmp
