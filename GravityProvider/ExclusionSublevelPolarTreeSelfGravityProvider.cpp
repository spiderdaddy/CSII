//
// Created by lf on 23/03/18.
//
#include <cstdlib>
#include <iostream>
#include <cmath>

#include "ExclusionSublevelPolarTreeSelfGravityProvider.h"


bool ExclusionSublevelPolarTreeSelfGravityProvider::isExcluded(QTNode *cell, QTNode *excluded_node) {
    if (cell == excluded_node) {
        return true;
    }
    for (QTNode *e : excluded_node->neighbour) {
        if (e == cell) {
            return true;
        }
    }
    return false;
}

void ExclusionSublevelPolarTreeSelfGravityProvider::calculate() {

    startTimer();

    disk->calcTreeValues(nullptr, disk->getQuadTree()->getHead());
    disk->getQuadTree()->printNodes();

    std::vector<Disk::Segment> &newSegment = *GravityProvider::pNewSegment;
    std::vector<Disk::Segment> &segment = *GravityProvider::pSegment;

    // Calculate acceleration for all points
    for (int r = 0; r < GravityProvider::num_radial_cells; r++) {
        for (int t = 0; t < GravityProvider::num_azimuthal_cells; t++) {

            int point = disk->getCellIndex(r, t);

            double ar = 0;
            double at = 0;

            QTNode *child = segment[point].node;
            QTNode *node = child->parent;

            // For the lowest level calculate all the leaves in the same parent
            for (int i = 0; i < 4; i++) {
                if ((node->leaf[i] != nullptr) &&
                    (node->leaf[i] != child)) {

                    calcGravityForNode(point, node->leaf[i], segment, ar, at);

                }
            }

            // Calculate all lowest level in neighbours of the parent
            for (QTNode *cell : node->neighbour) {
                for (int i = 0; i < 4; i++) {
                    if (cell->leaf[i] != nullptr) {

                        calcGravityForNode(point, cell->leaf[i], segment, ar, at);

                    }
                }
            }

            // for levels up to the resolution, calculate all the leaves of the neighbours
            // but not the leaves which were neighbours of the child
            child = node;
            node = node->parent;

            while (node->tree_level < resolution) {

                for (QTNode *cell : node->neighbour) {
                    for (int i = 0; i < 4; i++) {
                        if ((cell->leaf[i] != nullptr) &&
                            !isExcluded(cell->leaf[i], child)) {

                            calcGravityForNodeDepth(point, cell->leaf[i], segment, ar, at, child);
                        }
                    }
                }
                child = node;
                node = node->parent;
            }

            // Now calc everything else at the level we are at
            for (QTNode *cell : disk->getQuadTree()->getLevelVector(node->tree_level)) {
                for (int i = 0; i < 4; i++) {
                    if ((cell->leaf[i] != nullptr) &&
                        !isExcluded(cell->leaf[i], child)) {

                        calcGravityForNodeDepth(point, cell->leaf[i], segment, ar, at, child);
                    }
                }
            }

            newSegment[point].ar += ar;
            newSegment[point].at += at;
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

void ExclusionSublevelPolarTreeSelfGravityProvider::calcGravityForNodeDepth(
        int point,
        QTNode *node,
        vector<Disk::Segment> &segment,
        double &ar,
        double &at,
        QTNode *excluded_node) {

    /* Currently we assume that depth is either 0 or 1
     * I am proposing that this is equivalent to calculating the 0th or 1st
     * order derivative of the cell.
     *
     * For deeper levels we would implement a recursion down the tree for 'depth' levels
     */
    if (depth > 0) {
        for (int j = 0; j < 4; j++) {
            if (node->leaf[j] != nullptr) {

                calcGravityForNode(point, node->leaf[j], segment, ar, at);
            }
        }
    } else {
        calcGravityForNode(point, node, segment, ar, at);
    }
}


void ExclusionSublevelPolarTreeSelfGravityProvider::calcGravityForNode(
        int point,
        QTNode *node,
        vector<Disk::Segment> &segment,
        double &ar,
        double &at) {
    double diff_theta = segment[point].theta - node->theta;
    double cos_diff_theta = cos(diff_theta);
    double diff_r = segment[point].r - node->r;
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

