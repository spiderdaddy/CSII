//
// Created by lf on 23/03/18.
//
#include <cstdlib>
#include <iostream>
#include <cmath>

#include "SimplePolarTreeSelfGravityProvider.h"


void SimplePolarTreeSelfGravityProvider::calculate() {

    startTimer();

    disk->calcTreeValues(nullptr, disk->getQuadTree()->getHead());
    //disk->getQuadTree()->printNodes();

    std::vector<Disk::Segment> &newSegment = *GravityProvider::pNewSegment;
    std::vector<Disk::Segment> &segment = *GravityProvider::pSegment;

    // Calculate acceleration for all points
    for (int r = 0; r < GravityProvider::num_radial_cells; r++) {
        for (int t = 0; t < GravityProvider::num_azimuthal_cells; t++) {

            int i1 = disk->getCellIndex(r, t);

            double ar = 0;
            double at = 0;

            QTNode *child = segment[i1].node;
            QTNode *node = child->parent;

            while (node != nullptr) {

                for (int i = 0; i < 4; i++) {
                    if ((node->leaf[i] != nullptr) &&
                        (node->leaf[i] != child)) {

                        calcGravityForNode(i1, node->leaf[i], segment, ar, at);

                    }
                }
                child = node;
                node = node->parent;
            }

            newSegment[i1].ar += ar;
            newSegment[i1].at += at;

        }
        fprintf(
                stdout,
                "INFO: r=%4d\r",
                r
        );
        fflush(stdout);
    }

    stopTimer();
}

void SimplePolarTreeSelfGravityProvider::calcGravityForNode(int i, QTNode *node, vector<Disk::Segment> &segment,
                                                      double &ar, double &at) const {
    double diff_theta = segment[i].theta - node->theta;
    double cos_diff_theta = cos(diff_theta);
    double diff_r = segment[i].r - node->r;
    double exp_diff_r = exp(diff_r);
    double distance = (exp_diff_r * exp_diff_r)
                      + 1
                      - (2 * exp_diff_r * cos_diff_theta);
    distance = distance * sqrt(abs(distance));

    double a = 0.0;
    if (distance != 0.0) {
        a = -G * node->density * node->area / distance;
    }
    ar += a * (exp_diff_r - cos_diff_theta);
    at += a * sin(diff_theta);
}