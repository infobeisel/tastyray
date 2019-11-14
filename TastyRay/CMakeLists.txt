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


if(NOT TARGET OpenMP::OpenMP_CXX)

    find_package(Threads REQUIRED)
    add_library(OpenMP::OpenMP_CXX IMPORTED INTERFACE)
    set_property(TARGET OpenMP::OpenMP_CXX
                 PROPERTY INTERFACE_COMPILE_OPTIONS ${OpenMP_CXX_FLAGS})
    # Only works if the same flag is passed to the linker; use CMake 3.9+ otherwise (Intel, AppleClang)
    set_property(TARGET OpenMP::OpenMP_CXX
                 PROPERTY INTERFACE_LINK_LIBRARIES ${OpenMP_CXX_FLAGS} Threads::Threads)

endif()



#library
#add_library(TastyQuad STATIC
#executable
add_executable(TastyRay 
main.cpp
)

target_link_libraries (TastyRay PUBLIC
TastyQuad
			${SQLite3_LIBRARY}
GL					
			GLEW
			glfw3
X11 pthread dl sqlite3 gomp
)


target_compile_options(TastyRay PUBLIC -O3 -fopenmp)


