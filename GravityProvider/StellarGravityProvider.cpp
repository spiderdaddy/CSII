//
// Created by lf on 23/03/18.
//

#include "StellarGravityProvider.h"

void StellarGravityProvider::calculate()
{
    std::vector<Segment>& newSegment = *pNewSegment;
    std::vector<Segment>& segment = *pSegment;

    // Calculate acceleration for all points
    for (std::size_t r = 0; r < NUM_RADIAL_CELLS; r++) {
        for (std::size_t t = 0; t < NUM_AZIMUTHAL_CELLS; t++) {

            std::size_t i1 = (t * NUM_RADIAL_CELLS) + r;

            newSegment[i1].ar += -1.0 * *stellar_mass * G / (segment[i1].r * segment[i1].r);
            // newSegment[i1].at += 0;
        }
    }
}