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
include_directories("${PROJECT_SOURCE_DIR}/../../tinygltf")
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

target_link_libraries (TastyRay 
			"${PROJECT_SOURCE_DIR}/../../tastyquad/build/libTastyQuad.a"
						${SQLite3_LIBRARY}
GL					
			GLEW
pthread gomp
			${GLFW_LIBRARIES}

)


target_compile_options(TastyRay PUBLIC -O3 -fopenmp)


