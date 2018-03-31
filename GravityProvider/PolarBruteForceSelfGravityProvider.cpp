//
// Created by lf on 23/03/18.
//
#include <cstdlib>
#include <iostream>
#include <cmath>

#include "PolarBruteForceSelfGravityProvider.h"


void PolarBruteForceSelfGravityProvider::calculate() {

    int nearestNeighbours = max(GravityProvider::num_radial_cells, GravityProvider::num_azimuthal_cells) + 1;
    std::vector<Disk::Segment> &newSegment = *GravityProvider::pNewSegment;
    std::vector<Disk::Segment> &segment = *GravityProvider::pSegment;

    // Calculate acceleration for all points
    for (int r = 0; r < GravityProvider::num_radial_cells; r++) {
        for (int t = 0; t < GravityProvider::num_azimuthal_cells; t++) {

            int i1 = (r * GravityProvider::num_azimuthal_cells) + t;

            double ar = 0;
            double at = 0;

            int min_r = std::max(0, r - nearestNeighbours);
            int max_r = min(r + nearestNeighbours, GravityProvider::num_radial_cells);
            int  min_t = std::max(0, t - nearestNeighbours);
            int  max_t = min(t + nearestNeighbours, GravityProvider::num_azimuthal_cells);
            for (int sr = min_r; sr < max_r; sr++) {
                for (int st = min_t; st < max_t; st++) {

                    int i2 = (sr * GravityProvider::num_azimuthal_cells) + st;

                    if (i2 != i1) {
                        double diff_theta = segment[i1].theta - segment[i2].theta;
                        double cos_diff_theta = cos(diff_theta);
                        double diff_r = segment[i1].r - segment[i2].r;
                        double exp_diff_r = exp(diff_r);
                        double distance = (exp_diff_r * exp_diff_r)
                                          + 1
                                          - (2 * exp_diff_r * cos_diff_theta);
                        distance = distance * sqrt(abs(distance));
                        fflush(stdout);

                        double a = 0.0;
                        if (distance != 0.0) {
                            a = -G * segment[i2].density * segment[i2].area / distance;
                        }
                        ar += a * (exp_diff_r - cos_diff_theta);
                        at += a * sin(diff_theta);
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
}