cmake_minimum_required(VERSION 3.9)
project(src)

set(CMAKE_CXX_STANDARD 17)

set(SOURCE_FILES main.cpp)
add_executable(ewaldis ${SOURCE_FILES})

include_directories(src)

file(COPY data DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_subdirectory(src)
#add_subdirectory(ewaldis_tests)

target_link_libraries(ewaldis src libstdc++fs.a)