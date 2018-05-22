//
// Created by lf on 23/03/18.
//

#ifndef CSII_EXCLUSIONPOLARTREESELFGRAVITYPROVIDER_H
#define CSII_EXCLUSIONPOLARTREESELFGRAVITYPROVIDER_H


#include "GravityProvider.h"

class ExclusionPolarTreeSelfGravityProvider : public GravityProvider {

public:

    ExclusionPolarTreeSelfGravityProvider (Disk * d ) : GravityProvider ( d ) {}

    void calculate();

    void calcGravityForNode(int i1, QTNode *node, vector<Disk::Segment> &segment, double &ar, double &at) const;
    void calcGravityLeaf(int i1, QTNode *node, vector<Disk::Segment> &segment, double &ar, double &at, int level) const;

private:
    void calcGravityForNode0(int i1, QTNode *node, vector<Disk::Segment> &segment, double &ar, double &at) const;
    void calcGravityForNode1(int i1, QTNode *node, vector<Disk::Segment> &segment, double &ar, double &at) const;
};


#endif //CSII_EXCLUSIONPOLARTREESELFGRAVITYPROVIDER_H
