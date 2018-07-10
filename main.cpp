
#include <cstdlib>
#include <thread>
#include <zconf.h>

#include "disk.h"
#include "graphics.h"


#include "QuadTree.h"

void Initialize(int argc, char *argv[]) {

    InitializeGraphics(argc, argv);

}


int main(int argc, char *argv[]) {

    Initialize(argc, argv);

    // std::thread graphics (GraphicsMainLoop);

    //
    // Vortex
    //

    ApplyGravities(
            "/data/UZH/CSII/data1",
            "density",
            128,
            256,
            0,
            0,
            0,
            0);

    ApplyGravities(
            "/data/UZH/CSII/data1",
            "density",
            128,
            256,
            0,
            2,
            4,
            7);

    //
    // Planet
    //

    ApplyGravities(
            "/data/UZH/CSII/data1",
            "density_planet",
            128,
            256,
            0,
            0,
            0,
            0);

    ApplyGravities(
            "/data/UZH/CSII/data1",
            "density_planet",
            128,
            256,
            0,
            2,
            4,
            7);

    //
    //
    //

    ApplyGravities(
            "/data/UZH/CSII/data1",
            "density_final_05k",
            512,
            512,
            0,
            0,
            0,
            0);

    ApplyGravities(
            "/data/UZH/CSII/data1",
            "density_final_05k",
            512,
            512,
            0,
            2,
            4,
            7);

    //
    //
    //

    ApplyGravities(
            "/data/UZH/CSII/data1",
            "density_final_1k",
            1024,
            1024,
            0,
            2,
            4,
            7);

    ApplyGravities(
            "/data/UZH/CSII/data1",
            "density_final_1k",
            1024,
            1024,
            0,
            0,
            0,
            0);

    exit(EXIT_SUCCESS);
}







