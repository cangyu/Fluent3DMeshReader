cmake_minimum_required(VERSION 3.10)

project(MeshReader VERSION 1.0.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(inc)

add_library(VEC src/vec.cc)

add_library(XF src/xf.cc)
target_link_libraries(XF PUBLIC VEC)

add_library(REP src/rep.cc)
target_link_libraries(REP PRIVATE XF)

add_executable(REPRODUCE app/reproduce.cc)
target_link_libraries(REPRODUCE PUBLIC XF)

add_executable(ASCII_TRANSLATOR app/ascii.cc)
target_link_libraries(ASCII_TRANSLATOR PUBLIC REP)

add_executable(BINARY_TRANSLATOR app/binary.cc)
target_link_libraries(BINARY_TRANSLATOR PUBLIC REP)

install(TARGETS REPRODUCE
        RUNTIME DESTINATION ${CMAKE_SOURCE_DIR}/bin)

install(TARGETS ASCII_TRANSLATOR
        RUNTIME DESTINATION ${CMAKE_SOURCE_DIR}/bin)

install(TARGETS BINARY_TRANSLATOR
        RUNTIME DESTINATION ${CMAKE_SOURCE_DIR}/bin)
