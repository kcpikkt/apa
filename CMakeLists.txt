cmake_minimum_required(VERSION 3.13.1)

project(bigint VERSION 1.0.0)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

add_executable(example ${CMAKE_SOURCE_DIR}/example/main.cpp )
target_compile_options(example PUBLIC -Wall -Wextra -Wpedantic -g)
target_include_directories(example PUBLIC ${CMAKE_SOURCE_DIR}/include/)

add_executable(test ${CMAKE_SOURCE_DIR}/test/test.cpp )
target_compile_options(test PUBLIC -Wall -Wextra -Wpedantic -Wkeyword-macro -g)
target_link_libraries(test gmp)
target_include_directories(test PUBLIC ${CMAKE_SOURCE_DIR}/include/)
