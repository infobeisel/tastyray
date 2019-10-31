cmake_minimum_required (VERSION 3.1)
project (TastyRay)
set (CMAKE_CXX_STANDARD 14)
set (OpenGL_GL_PREFERENCE "GLVND")
find_package(glfw3)
find_package(OpenGL REQUIRED)
#find_package(GLEW REQUIRED)
include_directories("/usr/include/GL")
#find_package (SQLite3 REQUIRED)
find_package(GLFW REQUIRED)
include_directories("${PROJECT_SOURCE_DIR}")
include_directories(${GLFW_INCLUDE_DIRS})
include_directories("~/include/")
link_directories("~/lib/")
include_directories("~/usr/include/")
link_directories("~/usr/lib/")
link_directories("~/usr/lib/x86_64-linux-gnu/")

include_directories("${PROJECT_SOURCE_DIR}/../../tastyquad/TastyQuad")

#library
#add_library(TastyQuad STATIC
#executable
add_executable(TastyRay 
main.cpp
)

target_link_libraries (TastyRay 
TastyQuad
			${SQLite3_LIBRARY}
GL					
			GLEW
			glfw3
X11 pthread dl sqlite3
)


target_compile_options(TastyRay PUBLIC -O3)


