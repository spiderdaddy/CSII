//
// Created by lf on 23/03/18.
//
#include <cstdlib>
#include <iostream>
#include <cmath>

#include "PolarTreeSelfGravityProvider.h"


void PolarTreeSelfGravityProvider::calculate() {

    startTimer();

    disk->calcTreeValues(nullptr,  disk->getQuadTree()->getHead() );
    disk->getQuadTree()->printNodes();

    std::vector<Disk::Segment> &newSegment = *GravityProvider::pNewSegment;
    std::vector<Disk::Segment> &segment = *GravityProvider::pSegment;

    // Calculate acceleration for all points
    for (int r = 0; r < GravityProvider::num_radial_cells; r++) {
        for (int t = 0; t < GravityProvider::num_azimuthal_cells; t++) {

            int i1 = disk->getCellIndex(r, t);

            double ar = 0;
            double at = 0;

            QTNode * child = segment[i1].node;
            QTNode * node = child->parent;

            //while( node != nullptr ) {
            for( int j = 0; j < 3; j++ ) {

                for (int i = 0; i<4; i++) {
                    if ( ( node->leaf[i] != nullptr ) &&
                            (node->leaf[i] != child ) ) {

                        double diff_theta = segment[i1].theta - node->leaf[i]->theta;
                        double cos_diff_theta = cos(diff_theta);
                        double diff_r = segment[i1].r - node->leaf[i]->r;
                        double exp_diff_r = exp(diff_r);
                        double distance = (exp_diff_r * exp_diff_r)
                                          + 1
                                          - (2 * exp_diff_r * cos_diff_theta);
                        distance = distance * sqrt(abs(distance));

                        double a = 0.0;
                        if (distance != 0.0) {
                            a = -G * node->leaf[i]->density * node->leaf[i]->area / distance;
                        }
                        ar += a * (exp_diff_r - cos_diff_theta);
                        at += a * sin(diff_theta);


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