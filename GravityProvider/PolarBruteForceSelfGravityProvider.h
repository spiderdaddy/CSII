//
// Created by lf on 23/03/18.
//

#ifndef CSII_POLARBRUTEFORCESELFGRAVITYPROVIDER_H
#define CSII_POLARBRUTEFORCESELFGRAVITYPROVIDER_H


#include "GravityProvider.h"

class PolarBruteForceSelfGravityProvider : public GravityProvider {

public:

    PolarBruteForceSelfGravityProvider ( Disk * d ) : GravityProvider ( d ) {}

    void calculate();
};


#endif //CSII_POLARBRUTEFORCESELFGRAVITYPROVIDER_H
