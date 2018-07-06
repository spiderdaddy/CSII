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

            QTNode *node = segment[point].node;
            QTNode * child = nullptr;

            // For level 0 calculate all lowest level in neighbours of the node
            for (QTNode *cell : node->neighbour) {

                calcGravityForNode(point, cell, segment, ar, at);

            }

            child = node;
            node = node->parent;

            while (node->tree_level < resolution) {

                // Calc all cells at one level lower
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

/*void sum_mass( QTNode * root, QTNode * node, double &X, double &t ) {
    if ( node->tree_level == 0 ) {
        X += (node->r - root->r) * node->density * node->area;
        t += (node->theta - root->theta) * node->density * node->area;
    } else {
        for (int child = 0; child < 4; child++) {
            if (node->leaf[child] != nullptr) {
                sum_mass(root, node->leaf[child], X, t);
            }
        }
    }
}
*/
void ExclusionDifferentialPolarTreeSelfGravityProvider::calcGravityForNode(
        int point,
        QTNode *node,
        vector<Disk::Segment> &segment,
        double &ar,
        double &at) {

    // Calculate integrals for X_tilde and theta_tilde
    double X_mass = 0;
    double theta_mass = 0;
    if (depth >= 1) {

        double m;
        for (int child = 0; child < 4; child++) {
            if (node->leaf[child] != nullptr) {
                m = node->leaf[child]->density * node->leaf[child]->area;
                X_mass += (node->leaf[child]->r - node->r) * m;
                theta_mass += (node->leaf[child]->theta - node->theta) * m;
            }
        }

    }
    double X2_mass = 0;
    double theta2_mass = 0;
    double Xtheta_mass = 0;
    if (depth >= 2) {

        double m;
        for (int child = 0; child < 4; child++) {
            if (node->leaf[child] != nullptr) {
                m = node->leaf[child]->density * node->leaf[child]->area;
                X2_mass += (node->leaf[child]->r - node->r) * (node->leaf[child]->r - node->r) * m;
                theta2_mass += (node->leaf[child]->theta - node->theta) * (node->leaf[child]->theta - node->theta) * m;
                Xtheta_mass += (node->leaf[child]->r - node->r) * (node->leaf[child]->theta - node->theta) * m;
            }
        }

    }


    // The integral sigma * dXdT
    double mass = node->density * node->area;

    double delta_theta = segment[point].theta - node->theta;

    double delta_X = segment[point].r - node->r;

    double exp_delta_X = exp(delta_X);
    double cos_delta_theta = cos(delta_theta);
    double sin_delta_theta = sin(delta_theta);

    double l = 1.0 + (exp_delta_X * (exp_delta_X - (2.0 * cos_delta_theta)));
    double l32 = l * sqrt(abs(l));
    double l52 = l * l32;
    double l72 = l * l52;

    double dlX = 3.0 * exp_delta_X * (exp_delta_X - cos_delta_theta) / l52;
    double dltheta = 3.0 * exp_delta_X * sin_delta_theta / l52;

    double fr0 = (exp_delta_X - cos_delta_theta);
    double ftheta0 = sin_delta_theta;

    //
    // 0th order component
    //
    ar -= fr0 / l32 * mass;

    //
    // 1st order components
    //
    if (depth >= 1) {
        // X~ component of fr
        ar -= ((dlX * fr0) - (exp_delta_X / l32)) * X_mass;

        // theta~ component of fr
        ar -= ((dltheta * fr0) - (sin_delta_theta / l32)) * theta_mass;
    }


    //
    // 2nd order components
    //
    if (depth >= 2) {
        ar -= exp_delta_X * (
                (fr0 * 15.0 *
                 (exp_delta_X * exp_delta_X * (cos_delta_theta - exp_delta_X) * (cos_delta_theta - exp_delta_X))) / l72
                - (fr0 * 3.0 * (2.0 * exp_delta_X - cos_delta_theta)) / l52
                - (6.0 * exp_delta_X * exp_delta_X - exp_delta_X * cos_delta_theta) / l52
                + 1.0 / l32) * 0.5 * X2_mass;

        ar -= ((fr0 * 15.0 * (exp_delta_X * sin_delta_theta) * (exp_delta_X * sin_delta_theta)) / l72
               - (fr0 * 3.0 * exp_delta_X * cos_delta_theta) / l52
               - 6.0 * exp_delta_X * sin_delta_theta * sin_delta_theta / l52
               + cos_delta_theta / l32) * 0.5 * theta2_mass;

        ar -= (15.0 * (sin_delta_theta * (exp_delta_X * exp_delta_X * cos_delta_theta * cos_delta_theta +
                                        exp_delta_X * exp_delta_X * exp_delta_X * exp_delta_X)
                     - exp_delta_X * exp_delta_X * exp_delta_X * sin(2.0 * delta_theta)) / l72
               + 3.0 * (exp_delta_X * sin(2 * delta_theta) - 3.0 * exp_delta_X * exp_delta_X * sin_delta_theta) / l52) * Xtheta_mass;

    }

    //
    // 0th order component
    //
    at -= ftheta0 / l32 * mass;

    //
    // 1st order components
    //
    if (depth >= 1) {
        // X~ component of fr
        at -= dlX * ftheta0 * X_mass;

        // theta~ component of ftheta
        at -= ((dltheta * ftheta0) - (cos_delta_theta / l32)) * theta_mass;
    }


    //
    // 2nd order components
    //
    if (depth >= 2) {

        at -= (
                (15.0 * (exp_delta_X*(cos_delta_theta-exp_delta_X)) * (exp_delta_X*(cos_delta_theta-exp_delta_X)) * ftheta0 / l72 )
                -( 3.0 * exp_delta_X * ( 2.0 * exp_delta_X - cos_delta_theta) * ftheta0 / l52 )
              ) * 0.5 * X2_mass;

        at -= (
                      (15.0 * exp_delta_X*exp_delta_X*sin_delta_theta*sin_delta_theta*sin_delta_theta / l72 )
                - ( 3.0 * exp_delta_X * cos_delta_theta * sin_delta_theta / l52)
                -( 6.0 * exp_delta_X * cos_delta_theta * sin_delta_theta / l52 )
                -( sin_delta_theta / l32 )
              ) *0.5 * theta2_mass;

        at -= (
                      ( 1.5 * cos_delta_theta * ( -1.0* exp_delta_X * ( exp_delta_X - 2.0*cos_delta_theta - exp_delta_X*exp_delta_X)) / l52 )
        -( 7.5 * exp_delta_X*sin_delta_theta*sin_delta_theta * (-1.0 * exp_delta_X * ( exp_delta_X - 2.0*cos_delta_theta ) - exp_delta_X*exp_delta_X ) / l72 )
        - ( 3.0 * exp_delta_X * sin_delta_theta * sin_delta_theta / l52 )
              ) * Xtheta_mass;

    }

}


