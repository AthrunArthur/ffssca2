PROJECT(ffssca)
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
MESSAGE(STATUS "This is BINARY dir" ${PROJECT_BINARY_DIR})
MESSAGE(STATUS "This is SOURCE dir" ${PROJECT_SOURCE_DIR})
add_definitions(-std=c++11)

set(src bc_seq.cpp
  defs.h
  main.cpp
  read_data.cpp
  utils.h
util.cpp
)
include_directories(${PROJECT_SOURCE_DIR}/sprng2.0/include)
link_directories(${PROJECT_SOURCE_DIR}/sprng2.0/lib)
set(sprng_lib sprng m)

add_executable(${PROJECT_NAME} ${src})
target_link_libraries(${PROJECT_NAME} ${sprng_lib})
