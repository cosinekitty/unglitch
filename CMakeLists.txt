cmake_minimum_required (VERSION 3.0 FATAL_ERROR)
project(unglitch)

# Set high warning level and treat warnings as errors.
if(MSVC)
#    if(CMAKE_C_FLAGS MATCHES "/W[0-4]")
#        string(REGEX REPLACE "/W[0-4]" "/W4" CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
#    else()
#        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /W4")
#    endif()
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /W4 /WX")
elseif(CMAKE_COMPILER_IS_GNUCC)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wextra -Wall -Werror -pedantic -pedantic-errors -Wno-variadic-macros")
endif()

add_subdirectory(src)
