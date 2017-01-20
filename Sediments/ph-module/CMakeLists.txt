cmake_minimum_required (VERSION 2.8)
project (Reaction_parsing_module)

include_directories(/Users/MarkelovIgor/eigen)
include_directories(/Applications/MATLAB_R2014b.app/extern/include)
# include_directories(/usr/bin/include)
# find_package(Boost COMPONENTS unit_test_framework REQUIRED)
# set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

# find_package(MATLAB REQUIRED)


set(CMAKE_CXX_COMPILER "clang++") #compiler g++
# add_definitions("-std=c++11 -g -O0 -fopenmp -fPIC -Wall -Wextra -Wno-missing-field-initializers  -fdiagnostics-color=always -pthread")
add_definitions("-std=c++11 -O3 -ffast-math -ffinite-math-only -fopenmp -fPIC -Wall -Wextra -pthread -Wno-missing-field-initializers  -fdiagnostics-color=always")


add_library(imp_phreeqc STATIC IMPORTED)
set_property(TARGET imp_phreeqc PROPERTY
             IMPORTED_LOCATION /usr/local/lib/libiphreeqc.a)

add_executable(ph pH.cpp)
add_executable(pH_phreeqc pH_phreeqc.cpp)
add_executable(pH_num_dif pH_num_dif.cpp)

target_link_libraries(ph imp_phreeqc)
target_link_libraries(pH_phreeqc imp_phreeqc)
