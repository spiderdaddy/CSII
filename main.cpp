
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

    ApplyGravities();



    exit(EXIT_SUCCESS);
}







