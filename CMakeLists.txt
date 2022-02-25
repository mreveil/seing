
cmake_minimum_required(VERSION 3.4)

set(CMAKE_VERBOSE_MAKEFILE ON)

project(pyseing)

include_directories("${CMAKE_SOURCE_DIR}/src/")
# include_directories("${CMAKE_SOURCE_DIR}/python")

file(GLOB SOURCE_FILES "src/*")
# file (GLOB HEADER_FILES "src/*.h")

# Set up such that XCode organizes the files
# source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} FILES ${SOURCE_FILES} ${HEADER_FILES} ${PYTHON_FILES} )

# include(pybind11.cmake)
find_package(pybind11 REQUIRED)
pybind11_add_module(_pyseing
	${SOURCE_FILES}
)
