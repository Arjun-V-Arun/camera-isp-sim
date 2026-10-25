cmake_minimum_required(VERSION 3.16)
project(camera_isp_sim LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Default to an optimised build for single-config generators (Make, Ninja).
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

find_package(Threads REQUIRED)

add_library(isp STATIC
  src/image.cpp
  src/raw_meta.cpp
  src/stages.cpp
  src/pipeline.cpp
  src/analysis.cpp
  src/synth.cpp
)
target_include_directories(isp PUBLIC src)
target_link_libraries(isp PUBLIC Threads::Threads)

if(MSVC)
  target_compile_options(isp PUBLIC /W4 /permissive- /EHsc)
  target_compile_definitions(isp PUBLIC _CRT_SECURE_NO_WARNINGS NOMINMAX)
else()
  target_compile_options(isp PUBLIC -Wall -Wextra -Wpedantic)
endif()

add_executable(isp_sim src/main.cpp)
target_link_libraries(isp_sim PRIVATE isp)

enable_testing()
add_executable(isp_tests tests/test_main.cpp)
target_link_libraries(isp_tests PRIVATE isp)
add_test(NAME isp_tests COMMAND isp_tests)
