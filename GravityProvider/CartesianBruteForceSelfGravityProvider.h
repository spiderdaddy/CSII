//
// Created by lf on 23/03/18.
//

#ifndef CSII_CARTESIANBRUTEFORCESELFGRAVITYPROVIDER_H
#define CSII_CARTESIANBRUTEFORCESELFGRAVITYPROVIDER_H

#include "GravityProvider.h"

class CartesianBruteForceSelfGravityProvider : public GravityProvider {
public:

    CartesianBruteForceSelfGravityProvider ( Disk * d ) : GravityProvider ( d ) {}

    void calculate();

};

#endif //CSII_CARTESIANBRUTEFORCESELFGRAVITYPROVIDER_H
