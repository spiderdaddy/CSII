//
// Created by lf on 21/03/18.
//

#ifndef CSII_GRAVITYPROVIDER_H
#define CSII_GRAVITYPROVIDER_H

#include "../disk.h"
#include "../gravity.h"

class GravityProvider {

public:

    GravityProvider()= default;

    GravityProvider(
        std::vector<Segment> *pNewSegment,
        std::vector<Segment> *pSegment,
        int num_radial_cells,
        int num_azimuthal_cells,
        double * stellar_mass
    );

    virtual void calculate();

protected:
    std::vector<Segment> *pNewSegment;
    std::vector<Segment> *pSegment;
    int num_radial_cells;
    int num_azimuthal_cells;
    double * stellar_mass;
};


#endif //CSII_GRAVITYPROVIDER_H
