cmake_minimum_required (VERSION 3.1)
project (TastyRay)
set (CMAKE_CXX_STANDARD 17)
set (OpenGL_GL_PREFERENCE "GLVND")
find_package(PkgConfig REQUIRED)
pkg_search_module(GLFW REQUIRED glfw3)
find_package(OpenGL REQUIRED)
#find_package(GLEW REQUIRED)
include_directories("/usr/include/GL")
find_package (SQLite3 REQUIRED)
find_package(GLFW REQUIRED)
include_directories("${PROJECT_SOURCE_DIR}")
include_directories(${GLFW_INCLUDE_DIRS})

include_directories("${PROJECT_SOURCE_DIR}/../../tastyquad/TastyQuad")

#library
#add_library(TastyQuad STATIC
#executable
add_executable(TastyRay 
main.cpp
)

target_link_libraries (TastyRay 
			"${PROJECT_SOURCE_DIR}/../../tastyquad/ubuntuBuild/libTastyQuad.a"
						${SQLite3_LIBRARY}
GL					
			GLEW

			${GLFW_LIBRARIES}

)



