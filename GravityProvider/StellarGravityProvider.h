//
// Created by lf on 23/03/18.
//

#ifndef CSII_STELLARGRAVITYPROVIDER_H
#define CSII_STELLARGRAVITYPROVIDER_H

#include "GravityProvider.h"

class StellarGravityProvider : public GravityProvider {
public:

    StellarGravityProvider (
            std::vector<Segment> *pNS,
    std::vector<Segment> *pS,
    int num_r,
    int num_a,
    double * sm
    ) : GravityProvider ( pNS, pS, num_r, num_a, sm ) {}

    void calculate();
};

#endif //CSII_STELLARGRAVITYPROVIDER_H
