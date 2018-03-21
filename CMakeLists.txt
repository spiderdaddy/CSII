cmake_minimum_required(VERSION 3.9)
project(CSII)
set(CMAKE_CXX_STANDARD 11)

include_directories( BEFORE common )

add_executable(csii
        main.cpp
        graphics.cpp
        gravity.cpp
        disk.cpp
        #common/controls.cpp
        common/objloader.cpp
        common/quaternion_utils.cpp
        common/shader.cpp
        common/tangentspace.cpp
        #common/text2D.cpp
        #common/texture.cpp
        #common/vboindexer.cpp
        GravityProvider.cpp GravityProvider.h)

target_link_libraries(csii glut GL GLEW GLU )