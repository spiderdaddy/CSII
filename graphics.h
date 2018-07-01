//
// Created by lf on 20/03/18.
//

#ifndef CSII_GRAPHICS_H
#define CSII_GRAPHICS_H

#include "disk.h"
#include "gravity.h"

void RenderFunction(void);

void IdleFunction(void);

void TimerFunction(int);

void Cleanup(void);

void CreateVBO(void);

void DestroyVBO(void);

void CreateShaders(void);

void DestroyShaders(void);

void InitWindow(int argc, char *argv[]);

void ResizeFunction(int, int);

void InitializeGraphics(int argc, char *argv[]);

void setDisk(Disk *disk);

bool getKeyPressed();

void GraphicsMainLoop();

void saveImage( string filename );

void pause_render();
void resume_render();

bool isRendering();

#endif //CSII_GRAPHICS_H
