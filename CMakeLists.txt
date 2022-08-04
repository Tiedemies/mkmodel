cmake_minimum_required(VERSION 3.0.0)
project(MkModel VERSION 0.1.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -O3 -pedantic -std=gnu++17 -fopenmp")
SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fopenmp -lboost_date_time -std=gnu++17")

#Add library
add_subdirectory(utils)
add_subdirectory(simulators) 
add_subdirectory(mkmodel)
add_subdirectory(graphmodel)
add_subdirectory(algorithms)
add_subdirectory(data_gen)

#add executables
add_executable(RunEmpirical RunEmpirical.cpp)
add_executable(test_io test_io.cpp)
add_executable(vaccinator vaccinator.cpp)

install(TARGETS RunEmpirical DESTINATION build)
install(TARGETS test_io DESTINATION build)
install(TARGETS vaccinator DESTINATION build)

# Link libraries
target_link_libraries(RunEmpirical
    h_random
    defs
    ioroutines
    mkmodel    
    graphmodel
    mcsim
    optimizer
)

target_link_libraries(test_io
    stat_tester
    ioroutines
)

target_link_libraries(vaccinator
    defs
    ioroutines
    graphmodel
    hidden_cascade
    industry_cascade
)
