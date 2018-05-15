
#include <cstdlib>

#include "disk.h"
#include "graphics.h"

#include "QuadTree.h"

void Initialize(int argc, char *argv[]) {

    //Disk *disk = new Disk(128, 256, "/data/UZH/CSII/data1/density.data");
    Disk *disk = new Disk(128, 256, "/data/UZH/CSII/data1/density_planet.data");

    InitializeGraphics(argc, argv, disk);

}


int main(int argc, char *argv[]) {

//    Initialize(argc, argv);

//    GraphicsMainLoop();

    Disk *disk;

    disk = new Disk(128, 256, "/data/UZH/CSII/data1/density.data");
    ApplyGravity(disk, 4, "density1");
/*
    // do one to remove initial time
    ApplyGravity(disk, 7, "density0");
    ApplyGravity(disk, 0, "density1");
    ApplyGravity(disk, 0, "density2");
    ApplyGravity(disk, 0, "density3");
    ApplyGravity(disk, 0, "density4");

    for (int i = 3; i < disk->getQuadTree()->getMaxLevel() ; i++) {
        ApplyGravity(disk, i, "density");
    }

    disk = new Disk(128, 256, "/data/UZH/CSII/data1/density_planet.data");
    ApplyGravity(disk, 7, "density_planet0");
    ApplyGravity(disk, 0, "density_planet1");
    ApplyGravity(disk, 0, "density_planet2");
    ApplyGravity(disk, 0, "density_planet3");
    ApplyGravity(disk, 0, "density_planet4");

    for (int i = 3; i < disk->getQuadTree()->getMaxLevel() ; i++) {
        ApplyGravity(disk, i, "density_planet");
    }
*/



    exit(EXIT_SUCCESS);
}







