
#include <cstdlib>

#include "graphics.h"
#include "disk.h"


void Initialize(int argc, char *argv[]) {

    InitializeDisk();

    InitializeGraphics(argc, argv);

}


int main(int argc, char *argv[]) {

    Initialize(argc, argv);

    GraphicsMainLoop();

    exit(EXIT_SUCCESS);
}







