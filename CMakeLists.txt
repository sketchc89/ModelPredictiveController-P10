project(MPC)

cmake_minimum_required (VERSION 3.5)

# -g allows for gdb debugging
# turn on -03 for best performance
add_definitions(-std=c++11 -O3)

set(CXX_FLAGS)
# "-Wall")
set(CMAKE_CXX_FLAGS "${CXX_FLAGS}")

set(sources 
    src/MPC.cpp 
    src/main.cpp
    src/FG_eval.cpp)
set (tests 
    src/MPC_test.cpp
    src/FG_eval_test.cpp
    src/main_test.cpp)

include_directories(/usr/local/include)
link_directories(/usr/local/lib)
include_directories(src/Eigen-3.3)

add_executable(mpc ${sources})
target_link_libraries(mpc ipopt z ssl uv uWS)