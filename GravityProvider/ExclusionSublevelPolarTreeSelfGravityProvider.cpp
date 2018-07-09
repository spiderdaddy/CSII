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

            QTNode *node = segment[point].node;
            QTNode * child = nullptr;

            // For level 0 calculate all lowest level in neighbours of the node
            for (QTNode *cell : node->neighbour) {

                calcGravityForNode(point, cell, segment, ar, at);

            }

            child = node;
            node = node->parent;

            while (node->tree_level < resolution) {

                // Calc all children of the level we are at ignoring the exclusion zone of the previous level
                for (QTNode *cell : node->neighbour) {
                    for (int i = 0; i < 4; i++) {
                        if ((cell->leaf[i] != nullptr) &&
                            !isExcluded(cell->leaf[i], child)) {

                            calcGravityForNodeDepth(point, cell->leaf[i], segment, ar, at);
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

                        calcGravityForNodeDepth(point, cell->leaf[i], segment, ar, at);
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
        double &at) {

    /* Currently we assume that depth is either 0 or 1
     * I am proposing that this is equivalent to calculating the 0th or 1st
     * order derivative of the cell.
     *
     * For deeper levels we would implement a recursion down the tree for 'depth' levels
     */
    if (depth >= 1 && node->tree_level > 0 ) {
        for (int l1 = 0; l1 < 4; l1++) {
            if (node->leaf[l1] != nullptr) {
                if (depth >= 2 && node->tree_level > 1 ) {
                    for (int l2 = 0; l2 < 4; l2++) {
                        if (node->leaf[l1]->leaf[l2] != nullptr) {

                            calcGravityForNode(point, node->leaf[l1]->leaf[l2], segment, ar, at);
                        }
                    }
                } else {
                    calcGravityForNode(point, node->leaf[l1], segment, ar, at);
                }
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
    double delta_theta = segment[point].theta - node->theta;
    double cos_delta_theta = cos(delta_theta);
    double delta_r = segment[point].r - node->r;
    double exp_delta_r = exp(delta_r);
    double l = 1 + ( exp_delta_r * (exp_delta_r - (2 * cos_delta_theta)) );

    l = l * sqrt(abs(l));

    double mass = node->density * node->area;

    ar -= (exp_delta_r - cos_delta_theta) / l * mass;
    at -= sin(delta_theta) / l * mass;
}

