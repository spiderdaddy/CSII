//
// Created by lf on 23/03/18.
//

#ifndef CSII_POLARTREESELFGRAVITYPROVIDER_H
#define CSII_POLARTREESELFGRAVITYPROVIDER_H


#include "GravityProvider.h"

class SimplePolarTreeSelfGravityProvider : public GravityProvider {

public:

    SimplePolarTreeSelfGravityProvider (Disk * d ) : GravityProvider ( d ) {}

    void calculate();

    void calcGravityForNode(int i1, QTNode *node, vector<Disk::Segment> &segment, double &ar, double &at) const;
    void calcGravityLeaf(int i1, QTNode *node, vector<Disk::Segment> &segment, double &ar, double &at, int level) const;
};


#endif //CSII_POLARTREESELFGRAVITYPROVIDER_H
