//
// Created by lf on 23/03/18.
//
#include <cstdlib>
#include <iostream>
#include <cmath>

#include "ExclusionPolarTreeSelfGravityProvider.h"


bool isExcluded( QTNode * cell, QTNode * node ) {
    if (cell == node) {
        return true;
    }
    for (QTNode *e : node->neighbour ) {
        if ( e == cell ) {
            return true;
        }
    }
    return false;
}

void ExclusionPolarTreeSelfGravityProvider::calculate() {

    startTimer();

    disk->calcTreeValues(nullptr, disk->getQuadTree()->getHead());
    disk->getQuadTree()->printNodes();

    std::vector<Disk::Segment> &newSegment = *GravityProvider::pNewSegment;
    std::vector<Disk::Segment> &segment = *GravityProvider::pSegment;

    // Calculate acceleration for all points
//    for (int r = GravityProvider::num_radial_cells / 2 - 10; r < GravityProvider::num_radial_cells / 2 +10; r++) {
//        for (int t = GravityProvider::num_azimuthal_cells * 7/8 - 10; t < GravityProvider::num_azimuthal_cells * 7/8 +10; t++) {
  for (int r = 0; r < GravityProvider::num_radial_cells; r++) {
      for (int t = 0; t < GravityProvider::num_azimuthal_cells; t++) {

            int i1 = disk->getCellIndex(r, t);

            double ar = 0;
            double at = 0;

            QTNode *child = segment[i1].node;
            QTNode *node = child->parent;

            // For the lowest level calculate all the leaves in the same parent
            for (int i = 0; i < 4; i++) {
                if ((node->leaf[i] != nullptr) &&
                    (node->leaf[i] != child)) {

                    calcGravityForNode(i1, node->leaf[i], segment, ar, at);

                }
            }

            // Calculate all lowest level in neighbours of the parent
            for (QTNode *cell : node->neighbour ) {
                for (int i = 0; i < 4; i++) {
                    if (cell->leaf[i] != nullptr) {

                        calcGravityForNode(i1, cell->leaf[i], segment, ar, at);

                    }
                }
            }

            // for levels up to the resolution, calculate all the leaves of the neighbours
            // but not the leaves which were neighbours of the child
            child = node;
            node = node->parent;

            int level = node->tree_level;

            while (level < resolution ) {

                // Calc all lowest level in neighbours of the parent
                for (QTNode *cell : node->neighbour ) {
                    for (int i = 0; i < 4; i++) {
                        if ( (cell->leaf[i] != nullptr) &&
                             !isExcluded(cell->leaf[i], child) )   {
                            for (int j = 0; j < 4; j++) {
                                if (cell->leaf[i]->leaf[j] != nullptr) {

                                    calcGravityForNode(i1, cell->leaf[i]->leaf[j], segment, ar, at);
                                }
                            }
                        }
                    }
                }


                child = node;
                node = node->parent;
                level = node->tree_level;
            }

            // Now calc everything else at the level we are at
            for (QTNode *cell : disk->getQuadTree()->getLevelVector(level) ) {
                for (int i = 0; i < 4; i++) {
                    if ( (cell->leaf[i] != nullptr) &&
                         !isExcluded(cell->leaf[i], child) )   {
                        for (int j = 0; j < 4; j++) {
                            if (cell->leaf[i]->leaf[j] != nullptr) {

                                calcGravityForNode(i1, cell->leaf[i]->leaf[j], segment, ar, at);
                            }
                        }
                    }
                }
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

void ExclusionPolarTreeSelfGravityProvider::calcGravityForNode(
        int i,
        QTNode *node,
        vector<Disk::Segment> &segment,
        double &ar,
        double &at) const {
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

