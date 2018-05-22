//
// Created by lf on 23/03/18.
//

#ifndef CSII_EXCLUSIONSUBLEVELPOLARTREESELFGRAVITYPROVIDER_H
#define CSII_EXCLUSIONSUBLEVELPOLARTREESELFGRAVITYPROVIDER_H


#include "GravityProvider.h"

class ExclusionSublevelPolarTreeSelfGravityProvider : public GravityProvider {

public:

    ExclusionSublevelPolarTreeSelfGravityProvider (Disk * d ) : GravityProvider ( d ) {}

    void calculate();


private:
    void calcGravityForNode(int i1, QTNode *node, vector<Disk::Segment> &segment, double &ar, double &at);
    void calcGravityForNodeDepth(int i1, QTNode *node, vector<Disk::Segment> &segment, double &ar, double &at, QTNode *excluded_node);
    bool isExcluded( QTNode * cell, QTNode * excluded_node );
};


#endif //CSII_EXCLUSIONSUBLEVELPOLARTREESELFGRAVITYPROVIDER_H
