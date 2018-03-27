
#include <cstdlib>

#include "graphics.h"
#include "disk.h"

#define NUM_RADIAL_CELLS 128
#define NUM_AZIMUTHAL_CELLS 256


void Initialize(int argc, char *argv[]) {

//    Disk *disk = new Disk(NUM_RADIAL_CELLS, NUM_AZIMUTHAL_CELLS, "/data/UZH/CSII/data1/density.data");
    Disk *disk = new Disk(NUM_RADIAL_CELLS, NUM_AZIMUTHAL_CELLS, "/data/UZH/CSII/data1/density_planet.data");

    InitializeGraphics(argc, argv, disk);

}


int main(int argc, char *argv[]) {

    Initialize(argc, argv);

    GraphicsMainLoop();

    exit(EXIT_SUCCESS);
}







