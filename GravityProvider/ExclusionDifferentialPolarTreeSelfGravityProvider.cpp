//
// Created by lf on 23/03/18.
//
#include <cstdlib>
#include <iostream>
#include <cmath>

#include "ExclusionDifferentialPolarTreeSelfGravityProvider.h"


bool ExclusionDifferentialPolarTreeSelfGravityProvider::isExcluded(QTNode *cell, QTNode *excluded_node) {
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

void ExclusionDifferentialPolarTreeSelfGravityProvider::calculate() {

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
                    (node->leaf[i] != child)) {  // Note, if point offset then this must be included?

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


            // for levels up to the resolution, calculate all the neighbours
            // but not the leaves which were neighbours of the child
            child = node;
            node = node->parent;

            while (node->tree_level < resolution) {

                // Calc all lowest level in neighbours of the parent
                for (QTNode *cell : node->neighbour) {
                    for (int i = 0; i < 4; i++) {
                        if ((cell->leaf[i] != nullptr) &&
                            !isExcluded(cell->leaf[i], child)) {

                            calcGravityForNode(point, cell->leaf[i], segment, ar, at);
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

                        calcGravityForNode(point, cell->leaf[i], segment, ar, at);
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

void ExclusionDifferentialPolarTreeSelfGravityProvider::calcGravityForNode(
        int point,
        QTNode *node,
        vector<Disk::Segment> &segment,
        double &ar,
        double &at) {

    double X_mass = 0;
    double theta_mass = 0;
    for (int child = 0; child < 4; child++) {
        if (node->leaf[child] != nullptr) {
            X_mass += (node->leaf[child]->r - node->r) * node->leaf[child]->density * node->leaf[child]->area;
            theta_mass +=
                    (node->leaf[child]->theta - node->theta) * node->leaf[child]->density * node->leaf[child]->area;

        }
    }

    // The integral sigma * dXdT
    double mass = node->density * node->area;

    double delta_theta = segment[point].theta - node->theta;

    double delta_X = segment[point].r - node->r;

    double exp_delta_X = exp(delta_X);
    double cos_delta_theta = cos(delta_theta);
    double sin_delta_theta = sin(delta_theta);

    double l = 1 + (exp_delta_X * (exp_delta_X - 2 * cos_delta_theta));
    double l32 = l * sqrt(abs(l));
    double l52 = l * l32;

    double dlX = 3 * exp_delta_X * (exp_delta_X - cos_delta_theta) / l52;
    double dltheta = 3 * exp_delta_X * sin_delta_theta / l52;

    double fr0 = (exp_delta_X - cos_delta_theta);
    double ftheta0 = sin_delta_theta;

    double lfr;

    //
    // 0th order component
    //
    lfr = fr0 / l32;
    ar -= lfr * mass;

    //
    // 1st order components
    //
    if( depth >= 1 ) {
        lfr = dlX * fr0;
        lfr -= exp_delta_X / l32;
        // Note we are assuming G=1
        ar -= lfr * X_mass;

        // theta~ component of fr
        lfr = dltheta * fr0;
        lfr -= sin_delta_theta / l32;
        ar -= lfr * theta_mass;
    }

    double lft;

    //
    // 0th order component
    //
    lft = ftheta0 / l32;
    at -= lft * mass;

    //
    // 1st order components
    //
    if( depth >= 1 ) {
        lft = dlX * ftheta0;
        at -= lft * X_mass;

        // theta~ component of ftheta
        lft = dltheta * ftheta0;
        lft -= cos_delta_theta / l32;
        at -= lft * theta_mass;
    }
}


/*
 void ExclusionDifferentialPolarTreeSelfGravityProvider::calcGravityForNode1(
        int point,
        QTNode *node,
        vector<Disk::Segment> &segment,
        double &ar,
        double &at) {

    int child = 0;
    while (node->leaf[child] == nullptr && child < 4 ) {
        child ++;
    }
    if (child == 4 ) {
        return;
    }

    // Values we can reuse in both calculations
    double theta_tilde = abs(node->leaf[child]->theta - node->theta);
    double delta_theta = segment[point].theta - node->theta;
    double diff_theta = delta_theta - theta_tilde;

    double X_tilde = abs(node->leaf[child]->r - node->r);
    double delta_X = segment[point].r - node->r;
    double diff_X = delta_X - X_tilde;

    double exp_delta_X = exp(delta_X);
    double cos_delta_theta = cos(delta_theta);
    double sin_delta_theta = sin(delta_theta);
    double exp_diff_X = exp(diff_X);
    double cos_diff_theta = cos(diff_theta);
    double sin_diff_theta = sin(diff_theta);

    double l = 1 + (exp_diff_X * (exp_diff_X - 2 * cos_diff_theta));
    double l32 = l * sqrt(abs(l));
    double l52 = l * l;

    double lX = 3 * exp_delta_X * (exp_delta_X - cos_delta_theta) / l52;
    double ltheta = 3 * exp_delta_X * sin_delta_theta / l52;

    for (child = 0; child < 4; child++) {
        if (node->leaf[child] != nullptr) {

            // X~ component of fr
            double lfr = lX * (exp_diff_X - cos_diff_theta);
            lfr -= exp_delta_X / l32;
            double a = -G * lfr * X_tilde * (node->leaf[child]->density * node->leaf[child]->area);

            // theta~ component of fr
            lfr = ltheta * (exp_diff_X - cos_diff_theta);
            lfr -= sin_delta_theta / l32;
            a += -G * lfr * theta_tilde  * (node->leaf[child]->density * node->leaf[child]->area);

            ar += a;

            // X~ component of ftheta
            double lft = lX * sin_diff_theta;
            a = -G * lft * X_tilde * (node->leaf[child]->density * node->leaf[child]->area);

            // theta~ component of ftheta
            lft = ltheta * sin_diff_theta;
            lft -= cos_delta_theta / l32;
            a += -G * lft * theta_tilde  * (node->leaf[child]->density * node->leaf[child]->area);

            at += a;
        }
    }
}


 */

void ExclusionDifferentialPolarTreeSelfGravityProvider::calcGravityForNode2(
        int point,
        QTNode *node,
        vector<Disk::Segment> &segment,
        double &ar,
        double &at) {

    int child = 0;
    while (node->leaf[child] == nullptr && child < 4) {
        child++;
    }
    if (child == 4) {
        return;
    }

    // Values we can reuse in both calculations
    double theta_tilde = abs(node->leaf[child]->theta - node->theta);
    double delta_theta = segment[point].theta - node->theta;
    double diff_theta = delta_theta - theta_tilde;

    double X_tilde = abs(node->leaf[child]->r - node->r);
    double delta_X = segment[point].r - node->r;
    double diff_X = delta_X - X_tilde;

    double exp_delta_X = exp(delta_X);
    double cos_delta_theta = cos(delta_theta);
    double sin_delta_theta = sin(delta_theta);
    double exp_diff_X = exp(diff_X);
    double cos_diff_theta = cos(diff_theta);
    double sin_diff_theta = sin(diff_theta);

    double l = 1 + (exp_diff_X * (exp_diff_X - 2 * cos_diff_theta));
    double l32 = l * sqrt(abs(l));
    double l52 = l * l;

    double lX = 3 * exp_delta_X * (exp_delta_X - cos_delta_theta) / l52;
    double ltheta = 3 * exp_delta_X * sin_delta_theta / l52;

    for (child = 0; child < 4; child++) {
        if (node->leaf[child] != nullptr) {

            // X~ component of fr
            double lfr = lX * (exp_diff_X - cos_diff_theta);
            lfr -= exp_delta_X / l32;
            double a = -G * lfr * X_tilde * (node->leaf[child]->density * node->leaf[child]->area);

            // theta~ component of fr
            lfr = ltheta * (exp_diff_X - cos_diff_theta);
            lfr -= sin_delta_theta / l32;
            a += -G * lfr * theta_tilde * (node->leaf[child]->density * node->leaf[child]->area);

            ar += a;

            // X~ component of ftheta
            double lft = lX * sin_diff_theta;
            a = -G * lft * X_tilde * (node->leaf[child]->density * node->leaf[child]->area);

            // theta~ component of ftheta
            lft = ltheta * sin_diff_theta;
            lft -= cos_delta_theta / l32;
            a += -G * lft * theta_tilde * (node->leaf[child]->density * node->leaf[child]->area);

            at += a;
        }
    }
}
