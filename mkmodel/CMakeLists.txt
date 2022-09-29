#Cmakelist for mkmodels

add_library(mkmodel mkmodel.cpp)

target_link_libraries(mkmodel
    h_random
    defs
)