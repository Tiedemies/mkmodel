#Cmake list for data_gen libraries
add_library(stat_tester stat_tester.cpp)

target_link_libraries(stat_tester
    mkmodel
    ioroutines
    mcsim
    h_random
    defs
)