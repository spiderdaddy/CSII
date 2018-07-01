//
// Created by lf on 23/03/18.
//

#ifndef CSII_EXCLUSIONDIFFERENTIALPOLARTREESELFGRAVITYPROVIDER_H
#define CSII_EXCLUSIONDIFFERENTIALPOLARTREESELFGRAVITYPROVIDER_H


#include "GravityProvider.h"

class ExclusionDifferentialPolarTreeSelfGravityProvider : public GravityProvider {

public:

    ExclusionDifferentialPolarTreeSelfGravityProvider (Disk * d ) : GravityProvider ( d ) {}

    void calculate();


private:
    void calcGravityForNode (int point, QTNode *node, vector<Disk::Segment> &segment, double &ar, double &at);
    void calcGravityForNode0(int point, QTNode *node, vector<Disk::Segment> &segment, double &ar, double &at);
    void calcGravityForNode1(int point, QTNode *node, vector<Disk::Segment> &segment, double &ar, double &at);
    void calcGravityForNode2(int point, QTNode *node, vector<Disk::Segment> &segment, double &ar, double &at);
    bool isExcluded( QTNode * cell, QTNode * node );
};


#endif //CSII_EXCLUSIONDIFFERENTIALPOLARTREESELFGRAVITYPROVIDER_H
