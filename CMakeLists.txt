project(Fluent3DMeshReader)

cmake_minimum_required(VERSION 3.10)
set(CMAKE_CXX_STANDARD 17)

include_directories(inc)

add_library(VEC src/vec.cc)

add_library(XF src/xf.cc)
target_link_libraries(XF PUBLIC VEC)

add_library(REP src/rep.cc)
target_link_libraries(REP PRIVATE XF)

add_executable(CASE1 test/sp1.cc)
target_link_libraries(CASE1 PUBLIC XF)

add_executable(CASE2 test/sp2.cc)
target_link_libraries(CASE2 PUBLIC REP)
