cmake_minimum_required(VERSION 3.0.0)
project(PID VERSION 0.1.0)

include_directories(
  $ENV{ROOTSYS}/include
)

find_library(ROOTLIB Core $ENV{ROOTSYS})

add_definitions(
  -pthread -std=c++11 -m64 -fPIC
)

add_library(PID_Def SHARED PID_Def.C PID_Det.C)

link_directories(${ROOTLIB})
