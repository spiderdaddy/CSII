//
// Created by lf on 23/03/18.
//
#include <cstdlib>
#include <iostream>
#include <cmath>

#include "PolarBruteForceSelfGravityProvider.h"


void PolarBruteForceSelfGravityProvider::calculate() {

    std::vector<Disk::Segment>& newSegment = *GravityProvider::pNewSegment;
    std::vector<Disk::Segment>& segment = *GravityProvider::pSegment;

    // Calculate acceleration for all points
    for (std::size_t r = 0; r < GravityProvider::num_radial_cells; r++) {
        for (std::size_t t = 0; t < GravityProvider::num_azimuthal_cells; t++) {

            std::size_t i1 = (r * GravityProvider::num_azimuthal_cells) + t;

            double ar = 0;
            double at = 0;

            for (std::size_t sr = 0; sr < GravityProvider::num_radial_cells; sr++) {
                for (std::size_t st = 0; st < GravityProvider::num_azimuthal_cells; st++) {

                    std::size_t i2 = (sr * GravityProvider::num_azimuthal_cells) + st;

                    if (i2 != i1) {
                        double diff_theta = segment[i1].theta - segment[i2].theta;
                        double cos_diff_theta = cos(diff_theta);
                        double diff_r = segment[i1].r - segment[i2].r;
                        double exp_diff_r = exp(diff_r);
                        double distance = pow(
                                exp_diff_r * exp_diff_r
                                + 1
                                - ( 2 * exp_diff_r * cos_diff_theta )
                                , 1.5 )
                        ;
                        double a = 0;
                        if (distance > 0) {
                            a = G * segment[i2].density / distance;
                        }
                        ar += a * ( exp_diff_r - cos_diff_theta );
                        ar += a * sin(diff_theta) ;
                    }
                }
            }

            newSegment[i1].ar += ar;
            newSegment[i1].at += at;

        }
        fprintf(
                stdout,
                "INFO: r=%4lu\r",
                r
        );
        fflush(stdout);
    }
}