cmake_minimum_required(VERSION 2.8)

project(BalloonTest)

# OpenGL
# find_package(OpenGL REQUIRED)
# include_directories(${OPENGL_INCLUDE_DIR})

# Eigen
#set(CMAKE_MODULE_PATH  ${PROJECT_SOURCE_DIR}/cmake)
#find_package(Eigen REQUIRED)
#include_directories( ${EIGEN_INCLUDE_DIR} )

# GLUT
find_package(GLUT)
if(NOT GLUT_FOUND)
	set(GLUT_INCLUDE_DIRS ${CMAKE_CURRENT_LIST_DIR}/../GLUT)
	find_library(GLUT_LIBRARY glut32 PATHS ${CMAKE_CURRENT_LIST_DIR}/../GLUT )
	message(${GLUT_LIBRARY})
	set(GLUT_LIBRARIES ${GLUT_LIBRARY})
endif(NOT GLUT_FOUND)
include_directories(${GLUT_INCLUDE_DIRS})

include("CMakeLists.shared")

add_definitions(-DIGL_NO_OPENGL)

file(GLOB SOURCES
		${PROJECT_SOURCE_DIR}/src/*.cpp
		${PROJECT_SOURCE_DIR}/src/*.h
	)

LIST(APPEND SOURCES ${YIMAGE_SOURCES})
	
add_executable(${PROJECT_NAME} ${SOURCES} ${SHARED_SOURCES})
target_link_libraries(${PROJECT_NAME} ${SHARED_LIBRARIES} ${GLUT_LIBRARIES} ${PNG_LIBRARY} ${EMBREE_LIBRARIES})
