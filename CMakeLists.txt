cmake_minimum_required(VERSION 3.10)

project(MeshReader VERSION 1.0.1)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(inc)

add_library(VEC src/vec.cc)

add_library(XF src/xf.cc)
target_link_libraries(XF PUBLIC VEC)

add_library(REP src/rep.cc)
target_link_libraries(REP PRIVATE XF)

add_executable(XF_REPRODUCE app/reproduce.cc)
target_link_libraries(XF_REPRODUCE PUBLIC XF)

add_executable(XF_TRANSLATOR_ASCII app/ascii.cc)
target_link_libraries(XF_TRANSLATOR_ASCII PUBLIC REP)

add_executable(XF_TRANSLATOR_BINARY app/binary.cc)
target_link_libraries(XF_TRANSLATOR_BINARY PUBLIC REP)

install(TARGETS XF_REPRODUCE
        RUNTIME DESTINATION ${CMAKE_SOURCE_DIR}/bin)

install(TARGETS XF_TRANSLATOR_ASCII
        RUNTIME DESTINATION ${CMAKE_SOURCE_DIR}/bin)

install(TARGETS XF_TRANSLATOR_BINARY
        RUNTIME DESTINATION ${CMAKE_SOURCE_DIR}/bin)
