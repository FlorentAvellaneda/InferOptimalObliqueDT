cmake_minimum_required(VERSION 3.12)

set(projectName LearnODT)

message("cmake for ${projectName}")
project(${projectName})

set(CMAKE_CXX_STANDARD 20)

set(public_libs


)

set(private_libs
	fmt::fmt
	DiverLibs
	MaLib
	MesData
	svm
	z
	cadical
	EvalMaxSAT
	gmpxx
	gmp

)

FetchContent_Declare(
  CLI11
  GIT_REPOSITORY https://github.com/CLIUtils/CLI11.git
  GIT_SHALLOW TRUE 
  GIT_TAG v2.4.2)
FetchContent_MakeAvailable(CLI11)
list(APPEND private_libs CLI11::CLI11)

file(
	GLOB_RECURSE
	source_files
	src/*
)

add_executable (
    ${projectName}
    ${source_files}
)

target_include_directories(${projectName} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)

target_link_libraries(${projectName} PUBLIC ${public_libs} PRIVATE ${private_libs})

