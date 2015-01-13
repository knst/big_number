project(big_number)

set(SOURCE
    additional.cpp
    BN.cpp
    BNsign.cpp
    testing.cpp
    main.cpp
)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -g")

add_executable(big_number ${SOURCE})
