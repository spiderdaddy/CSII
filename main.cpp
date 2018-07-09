
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

    std::thread graphics (GraphicsMainLoop);

    ApplyGravities(
            "/data/UZH/CSII/data1",
            "density",
            128,
            256,
            0,
            2,
            6,
            6);


    ApplyGravities(
            "/data/UZH/CSII/data1",
            "density_planet",
            128,
            256,
            0,
            2,
            6,
            6);

/*
    ApplyGravities(
            "/data/UZH/CSII/data1",
            "density_final_05k",
            512,
            512,
            0,
            1,
            3,
            7);
/*
    ApplyGravities(
            "/data/UZH/CSII/data1",
            "density_final_1k",
            1024,
            1024,
            0,
            1,
            1,
            6);
*/
    exit(EXIT_SUCCESS);
}







