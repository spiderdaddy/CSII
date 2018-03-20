//
// Created by lf on 20/03/18.
//

#ifndef CSII_GRAVITY_H
#define CSII_GRAVITY_H

#include "disk.h"

void ApplyBruteForceGravity(
        std::vector<Segment> newSegment,
        std::vector<Segment> segment,
        double *stellar_mass,
        double *escape_mass
);

//void ApplyXNearestNeighbourGravity();
//void ApplyNearestNeighbourGravity();

#define G 6.6e-11
#define M_EARTH   1e22
#define M_SUN     2e30
#define M_JUPITER 1e27
#define STAR
#define SELF

#endif //CSII_GRAVITY_H
