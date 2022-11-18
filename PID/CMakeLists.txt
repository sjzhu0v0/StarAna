cmake_minimum_required(VERSION 3.0.0)
project(PID VERSION 0.1.0)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# include_directories(
# $ENV{ROOTSYS}/include
# ${PROJECT_SOURCE_DIR}
# )
find_package(ROOT 6.18 CONFIG REQUIRED)

add_library(PID_Def SHARED PID_Def.C PID_Det.C)
link_directories(${ROOT_LIBRARY_DIR})

add_definitions(
  -pthread -m64 -fPIC -shared
)

target_include_directories(PID_Def PUBLIC
  ${ROOT_INCLUDE_DIRS}
  ${PROJECT_SOURCE_DIR}
)

message("ROOT_INCLUDE_DIRS" ${ROOT_INCLUDE_DIRS})
message("PROJECT_SOURCE_DIR" ${PROJECT_SOURCE_DIR})

target_link_libraries(PID_Def PUBLIC
  ${ROOT_LIBRARIES} -L${ROOT_LIBRARY_DIR}
)