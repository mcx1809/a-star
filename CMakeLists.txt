cmake_minimum_required(VERSION 3.5)

project(test)

find_package(Boost 1.71)
add_library(boost INTERFACE)
target_include_directories(boost INTERFACE ${Boost_INCLUDE_DIRS})
target_link_libraries(boost INTERFACE ${Boost_LIBRARIES})

add_executable(${CMAKE_PROJECT_NAME} src/test.cc)
target_include_directories(${CMAKE_PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(${CMAKE_PROJECT_NAME} boost)
