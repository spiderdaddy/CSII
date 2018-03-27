//
// Created by lf on 23/03/18.
//

#ifndef CSII_CARTESIANBRUTEFORCESELFGRAVITYPROVIDER_H
#define CSII_CARTESIANBRUTEFORCESELFGRAVITYPROVIDER_H

#include "GravityProvider.h"

class CartesianBruteForceSelfGravityProvider : public GravityProvider {
public:

    CartesianBruteForceSelfGravityProvider (
            std::vector<Disk::Segment> *pNS,
    std::vector<Disk::Segment> *pS,
    int num_r,
    int num_a,
    double * sm
    ) : GravityProvider ( pNS, pS, num_r, num_a, sm ) {}

    void calculate();

};

#endif //CSII_CARTESIANBRUTEFORCESELFGRAVITYPROVIDER_H
