//
// Created by lf on 23/03/18.
//

#ifndef CSII_POLARTREESELFGRAVITYPROVIDER_H
#define CSII_POLARTREESELFGRAVITYPROVIDER_H


#include "GravityProvider.h"

class PolarTreeSelfGravityProvider : public GravityProvider {

public:

    PolarTreeSelfGravityProvider (Disk * d ) : GravityProvider ( d ) {}

    void calculate();
};


#endif //CSII_POLARTREESELFGRAVITYPROVIDER_H
