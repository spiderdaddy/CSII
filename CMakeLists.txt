cmake_minimum_required(VERSION 3.9)
project(CSII)
set(CMAKE_CXX_STANDARD 11)

include_directories( BEFORE common GravityProvider )

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
        GravityProvider/GravityProvider.cpp
        GravityProvider/PolarBruteForceSelfGravityProvider.cpp
        GravityProvider/SimplePolarTreeSelfGravityProvider.cpp
        GravityProvider/ExclusionSublevelPolarTreeSelfGravityProvider.cpp
        GravityProvider/ExclusionDifferentialPolarTreeSelfGravityProvider.cpp
        GravityProvider/StellarGravityProvider.cpp
        GravityProvider/CartesianBruteForceSelfGravityProvider.cpp

        QuadTree.cpp

        )

target_link_libraries(csii glut GL GLEW GLU )