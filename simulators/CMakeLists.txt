# Cmake list for directory Simulators
add_library(hidden_cascade hidden_cascade.cpp)
add_library(announcement_cascade announcement_cascade.cpp)
add_library(industry_cascade industry_cascade.cpp)
add_library(mcsim mcsim.cpp)

# Link libraries
target_link_libraries(hidden_cascade
    h_random
    defs
    ioroutines
    graphmodel   
)

target_link_libraries(industry_cascade
    h_random
    defs
    ioroutines
    graphmodel
    hidden_cascade   
)

target_link_libraries(mcsim
    graphmodel
    defs
    mkmodel
    h_random
)

target_link_libraries(announcement_cascade
    h_random
    defs
    ioroutines
    graphmodel
    hidden_cascade
    industry_cascade   
)

