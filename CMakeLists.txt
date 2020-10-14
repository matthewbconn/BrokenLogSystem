cmake_minimum_required(VERSION 3.17)
project(VivadoLogSystem)

set(CMAKE_CXX_STANDARD 14)
include_directories("C:/Xilinx/Vivado/2019.2/include")

add_executable(VivadoLogSystem main.cpp lognum.cpp lognum.h ConversionEngine.cpp)

