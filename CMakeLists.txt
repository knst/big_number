project(big_number)

set(SOURCE
    additional.cpp
    BN.cpp
    BNsign.cpp
    testing.cpp
    main.cpp
)

set(WARNING_FLAGS "-pedantic -Wall -Wextra -Wshadow -Wformat=2 -Wfloat-equal -Wlogical-op -Wcast-qual -Wcast-align")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${WARNING_FLAGS} -std=c++11 -g")
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D _GLIBCXX_DEBUG -fsanitize=address")

add_executable(big_number ${SOURCE})
