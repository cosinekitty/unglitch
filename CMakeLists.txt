cmake_minimum_required (VERSION 3.0 FATAL_ERROR)
project(unglitch)

# Set high warning level and treat warnings as errors.
if(MSVC)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /W4 /WX")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4 /WX")
elseif(CMAKE_COMPILER_IS_GNUCC)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wextra -Wall -Werror -pedantic -pedantic-errors -Wno-variadic-macros")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wextra -Wall -Werror -pedantic -pedantic-errors -Wno-variadic-macros")
endif()

add_subdirectory(src)
