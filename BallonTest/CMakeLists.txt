cmake_minimum_required(VERSION 2.8)

project(BalloonTest)

# OpenGL
find_package(OpenGL REQUIRED)
include_directories(${OPENGL_INCLUDE_DIR})

# GLUT
find_package(GLUT)
if(NOT GLUT_FOUND)
	set(GLUT_INCLUDE_DIRS ${CMAKE_CURRENT_LIST_DIR}/../GLUT)
	find_library(GLUT_LIBRARY glut32 PATHS ${CMAKE_CURRENT_LIST_DIR}/../GLUT )
	message(${GLUT_LIBRARY})
	set(GLUT_LIBRARIES ${GLUT_LIBRARY})
endif(NOT GLUT_FOUND)
include_directories(${GLUT_INCLUDE_DIRS})

# Eigen
set(CMAKE_MODULE_PATH  ${PROJECT_SOURCE_DIR}/cmake)
find_package(Eigen REQUIRED)
include_directories( ${EIGEN_INCLUDE_DIR} )

file(GLOB ex_files
		${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/src/*.h
	)
	
add_executable(BalloonTest ${ex_files})

target_link_libraries(BalloonTest ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES})
