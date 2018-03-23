//
// Created by lf on 21/03/18.
//

#include "GravityProvider.h"



GravityProvider::GravityProvider (
        std::vector<Segment> *pNS,
        std::vector<Segment> *pS,
        int num_r,
        int num_a,
        double * sm
) {
    pNewSegment = pNS;
    pSegment = pS;
    num_radial_cells = num_r;
    num_azimuthal_cells = num_a;
    stellar_mass = sm;
}

void GravityProvider::calculate() {

}



