//
// Created by lf on 20/03/18.
//

#ifndef CSII_GRAVITY_H
#define CSII_GRAVITY_H

#include "disk.h"
#include <GravityProvider.h>


void ApplyGravities();
void ApplyGravity(Disk *disk, string data_name, GravityProvider *gp, int resolution, int depth);

//void ApplyXNearestNeighbourGravity();
//void ApplyNearestNeighbourGravity();


#endif //CSII_GRAVITY_H
