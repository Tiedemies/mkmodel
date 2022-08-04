#CMakelist for algorithms
add_library(optimizer optimizer.cpp)

target_link_libraries(optimizer
    graphmodel
    mkmodel
    ioroutines
    defs
    h_random
    mcsim
)