//
// Created by lf on 23/03/18.
//

#ifndef CSII_STELLARGRAVITYPROVIDER_H
#define CSII_STELLARGRAVITYPROVIDER_H

#include "GravityProvider.h"

class StellarGravityProvider : public GravityProvider {
public:

    StellarGravityProvider ( Disk * d ) : GravityProvider ( d ) {}

    void calculate();
};

#endif //CSII_STELLARGRAVITYPROVIDER_H
