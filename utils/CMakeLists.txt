add_library(h_random h_random.cpp)
add_library(defs defs.cpp)
add_library(ioroutines ioroutines.cpp)

target_link_libraries(ioroutines
    defs
)