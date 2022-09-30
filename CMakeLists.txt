cmake_minimum_required(VERSION 3.20)
#project(MkModel VERSION 0.1.0)


set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

SET(CMAKE_CONFIGURATION_TYPES "Debug;Release" CACHE STRING "Configs" FORCE)
IF(DEFINED CMAKE_BUILD_TYPE AND CMAKE_VERSION VERSION_GREATER "3.0")
  SET_PROPERTY(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS  ${CMAKE_CONFIGURATION_TYPES})
ENDIF()

set(CMAKE_INSTALL_PREFIX .)

if(CMAKE_BUILD_TYPE MATCHES "Release")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -O3 -pedantic -std=gnu++17 -fopenmp -DNDEBUG")
  SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wall -O3 -fopenmp -lboost_date_time -std=gnu++17 -DNDEBUG")
endif()

if(CMAKE_BUILD_TYPE MATCHES "Debug")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -g3 -pedantic -std=gnu++17")
  SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lboost_date_time -std=gnu++17")
endif()


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

# Install into build-type directory
install(TARGETS RunEmpirical test_io vaccinator DESTINATION ${CMAKE_BUILD_TYPE})


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
    vacc_minimizer
)
