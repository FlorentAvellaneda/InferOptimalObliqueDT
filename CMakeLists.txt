cmake_minimum_required(VERSION 3.12)

set(projectName LearnODT)

message("cmake for ${projectName}")
project(
  ${projectName}
  VERSION 0.1
  DESCRIPTION "Projet ${projectName}"
  LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 20)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

include(FetchContent)

FetchContent_Declare(
  fmt
  GIT_REPOSITORY https://github.com/fmtlib/fmt.git
  GIT_SHALLOW TRUE 
  GIT_TAG 11.0.0)
FetchContent_MakeAvailable(fmt)
list(APPEND externLibs fmt::fmt)


add_subdirectory(lib/MaLib)
add_subdirectory(lib/cadical)
add_subdirectory(lib/EvalMaxSAT)
add_subdirectory(lib/MesData)
add_subdirectory(lib/DiverLibs)

add_subdirectory(main)



