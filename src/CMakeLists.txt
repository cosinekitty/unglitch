include_directories(
    ${CMAKE_SOURCE_DIR}/src
)

add_executable(unglitch
    ${CMAKE_SOURCE_DIR}/src/unglitch.cpp
    ${CMAKE_SOURCE_DIR}/src/tinyxml2.cpp
)

install(TARGETS unglitch DESTINATION ${CMAKE_SOURCE_DIR}/bin/${CMAKE_BUILD_TYPE})
