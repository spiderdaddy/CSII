//
// Created by lf on 23/03/18.
//

#include "StellarGravityProvider.h"

void StellarGravityProvider::calculate()
{
    std::vector<Disk::Segment>& newSegment = *pNewSegment;
    std::vector<Disk::Segment>& segment = *pSegment;

    // Calculate acceleration for all points
    for (std::size_t r = 0; r < num_radial_cells; r++) {
        for (std::size_t t = 0; t < num_azimuthal_cells; t++) {

            std::size_t i1 = (t * num_radial_cells) + r;

            newSegment[i1].ar += -1.0 * *stellar_mass * G / (segment[i1].r * segment[i1].r);
            // newSegment[i1].at += 0;
        }
    }
}