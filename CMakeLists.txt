project(big_number)
cmake_minimum_required(VERSION 2.8)

set(SOURCE
    additional.cpp
    BN.cpp
    BNsign.cpp
)

set(SOURCE_EXAMPLE
    testing.cpp
    main.cpp
)

set(WARNING_FLAGS "-pedantic -Wall -Wextra -Wshadow -Wformat=2 -Wfloat-equal -Wlogical-op -Wcast-qual -Wcast-align")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${WARNING_FLAGS} -std=c++11 -g")
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D _GLIBCXX_DEBUG -fsanitize=address")

add_library(big_number_static STATIC ${SOURCE})
add_library(big_number_shared SHARED ${SOURCE})
set_target_properties(big_number_static PROPERTIES OUTPUT big_number)
set_target_properties(big_number_shared PROPERTIES OUTPUT big_number)

add_executable(bignumber_example ${SOURCE_EXAMPLE})
target_link_libraries(bignumber_example big_number_static)
