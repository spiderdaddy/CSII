//
// Created by lf on 23/03/18.
//
#include <cstdlib>
#include <iostream>
#include <cmath>

#include "CartesianBruteForceSelfGravityProvider.h"

void CartesianBruteForceSelfGravityProvider::calculate() {

    std::vector<Segment>& newSegment = *pNewSegment;
    std::vector<Segment>& segment = *pSegment;

    // Calculate acceleration for all points
    for (std::size_t r = 0; r < NUM_RADIAL_CELLS; r++) {
        for (std::size_t t = 0; t < NUM_AZIMUTHAL_CELLS; t++) {

            std::size_t i1 = (r * NUM_AZIMUTHAL_CELLS) + t;

            double ar = 0;
            double at = 0;

            for (std::size_t sr = 0; sr < NUM_RADIAL_CELLS; sr++) {
                for (std::size_t st = 0; st < NUM_AZIMUTHAL_CELLS; st++) {

                    std::size_t i2 = (sr * NUM_AZIMUTHAL_CELLS) + st;

                    if (i2 != i1) {

                        double x = segment[i2].x - segment[i1].x;
                        double y = segment[i2].y - segment[i1].y;
                        double d = sqrt((x * x) * (y * y));
                        double alpha = atan(y/x);
                        double f = d * d;
                        double a = 0;
                        if (f > 0) {
                            a = G * segment[i2].density * segment[i2].area / f;
                        }
                        double dt = segment[i1].theta - alpha;
                        at += a * sin(dt) ;
                        ar += a * cos(dt) ;
                    }
                }
            }
            newSegment[i1].ar += ar;
            newSegment[i1].at += at;
        }
    }
}