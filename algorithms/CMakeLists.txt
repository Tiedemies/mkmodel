#CMakelist for algorithms
add_library(optimizer optimizer.cpp)
add_library(vacc_minimizer vacc_minimizer.cpp)

target_link_libraries(optimizer
    graphmodel
    mkmodel
    ioroutines
    defs
    h_random
    mcsim
)

target_link_libraries(vacc_minimizer
    graphmodel
    mkmodel
    ioroutines
    defs
    h_random
    industry_cascade
    hidden_cascade
    announcement_cascade
)