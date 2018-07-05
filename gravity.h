//
// Created by lf on 20/03/18.
//

#ifndef CSII_GRAVITY_H
#define CSII_GRAVITY_H

#include "disk.h"
#include <GravityProvider.h>


void ApplyGravities( string path,
                     string dataname,
                     unsigned r_cells,
                     unsigned theta_cells,
                     unsigned min_depth,
                     unsigned max_depth,
                     unsigned min_resolution,
                     unsigned max_resolution) ;

void ApplyGravity(Disk *disk, string pathname, string data_name, GravityProvider *gp, int resolution, int depth) ;

//void ApplyXNearestNeighbourGravity();
//void ApplyNearestNeighbourGravity();


#endif //CSII_GRAVITY_H
