cmake_minimum_required(VERSION 3.0.0)
project(PID VERSION 0.1.0)

include_directories(
  /data/Software/ROOT/6.24.00/include
)

find_library(ROOTLIB Core ${ROOTSYS})

add_definitions(
  -pthread -std=c++11 -m64
)

add_library(PID_Def PID_Def.C PID_Det.C)

link_directories(${ROOTLIB})

# target_link_directories(PID_Det
#   ${ROOTSYS}/lib
#   ./
# )

