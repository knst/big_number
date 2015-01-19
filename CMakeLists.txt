project(big_number)

set(SOURCE
    additional.cpp
    BN.cpp
    BNsign.cpp
    testing.cpp
    main.cpp
)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -g -Wall -Wextra -pedantic -Wshadow -Wformat=2 -Wfloat-equal -Wlogical-op -Wcast-qual -Wcast-align")

add_executable(big_number ${SOURCE})
