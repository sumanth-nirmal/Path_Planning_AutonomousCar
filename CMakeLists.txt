project(Path_Planning_Bosch)

cmake_minimum_required (VERSION 2.8.3)

add_definitions(-std=c++11)

set(CXX_FLAGS "-Wall")
set(CMAKE_CXX_FLAGS, "${CXX_FLAGS}")

include_directories(${CMAKE_SOURCE_DIR}/include)
set(sources src/main.cpp src/traj_planner.cpp)

if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin") 

include_directories(/usr/local/include)
include_directories(/usr/local/opt/openssl/include)
link_directories(/usr/local/lib)
link_directories(/usr/local/opt/openssl/lib)
link_directories(/usr/local/Cellar/libuv/1.11.0/lib)

endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin") 

add_executable(path_planning ${sources})

target_link_libraries(path_planning z ssl uv uWS)
