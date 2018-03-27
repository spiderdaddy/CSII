//
// Created by lf on 20/03/18.
//

#include <vector>
#include <cstdlib>
#include <algorithm>

#include <glm/glm.hpp>

#include "disk.h"
#include "gravity.h"
#include "CartesianBruteForceSelfGravityProvider.h"
#include "PolarBruteForceSelfGravityProvider.h"
#include "StellarGravityProvider.h"

double r2 = 1 / std::sqrt(2);

void initialiseAcceleration(
        std::vector<Disk::Segment> *pNewSegment,
        std::vector<Disk::Segment> *pSegment) {

    std::vector<Disk::Segment> &newSegment = *pNewSegment;
    std::vector<Disk::Segment> &segment = *pSegment;

    std::for_each(newSegment.begin(), newSegment.end(), [](Disk::Segment &el) {
        el.ar = 0;
        el.at = 0;
    });

    std::for_each(segment.begin(), segment.end(), [](Disk::Segment &el) {
        el.ar = 0;
        el.at = 0;
    });
}


void ApplyAcceleration(
        std::vector<Disk::Segment> *pNewSegment
) {

    std::vector<Disk::Segment> &newSegment = *pNewSegment;

    // apply acceleration to all velocities
    for (int i = 0; i < newSegment.size(); i++) {
        newSegment[i].vr += newSegment[i].ar;
        newSegment[i].vt += newSegment[i].at;
        newSegment[i].pr.clear();
        newSegment[i].pt.clear();
    }

}

void MoveMass( Disk *disk ) {
    std::vector<Disk::Segment> &newSegment = *disk->getNewSegment();

    // populate points
    double radius_step = (OUTER_RADIUS - INNER_RADIUS) / disk->get_num_radial_cells();
    double radius_step_2 = radius_step / 2;

    // based on velocities, move the matter, conserve momentum
    for (std::size_t r = 0; r < disk->get_num_radial_cells(); r++) {
        for (std::size_t t = 0; t < disk->get_num_azimuthal_cells(); t++) {

            std::size_t i1 = (t * disk->get_num_radial_cells()) + r;

            long *n = &newSegment[i1].n[0];

            double theta_step = 2 * (double) M_PI * newSegment[i1].r / disk->get_num_azimuthal_cells();
            double theta_step_2 = theta_step / 2;
            double r1 = newSegment[i1].r - radius_step_2;
            double r2 = r1 + newSegment[i1].vr;
            double t1 = newSegment[i1].theta - theta_step_2;
            double t2 = t1 + newSegment[i1].vt;

            double a1 = 0;
            double a2 = 0;
            double a3 = 0;
            int n1, n2, n3;
            if (r2 >= r1) {
                if (t2 >= t1) {
                    // Q1
                    a1 = glm::min(radius_step, r2 - r1) * glm::min(theta_step, t1 + theta_step - t2);
                    a2 = glm::min(radius_step, r2 - r1) * glm::min(theta_step, t2 - t1);
                    a3 = glm::min(radius_step, r1 + radius_step - r2) * glm::min(theta_step, t2 - t1);
                    n1 = 4;
                    n2 = 3;
                    n3 = 2;
                } else {
                    // Q4
                    a1 = glm::min(radius_step, r2 - r1) * glm::min(theta_step, t2 + theta_step - t1);
                    a2 = glm::min(radius_step, r2 - r1) * glm::min(theta_step, t1 - t2);
                    a3 = glm::min(radius_step, r1 + radius_step - r2) * glm::min(theta_step, t1 - t2);
                    n1 = 4;
                    n2 = 5;
                    n3 = 6;
                }
            } else {
                if (t2 >= t1) {
                    // Q2
                    a1 = glm::min(radius_step, r1 - r2) * glm::min(theta_step, t1 + theta_step - t2);
                    a2 = glm::min(radius_step, r1 - r2) * glm::min(theta_step, t2 - t1);
                    a3 = glm::min(radius_step, r2 + radius_step - r1) * glm::min(theta_step, t2 - t1);
                    n1 = 0;
                    n2 = 1;
                    n3 = 2;
                } else {
                    // Q3
                    a1 = glm::min(radius_step, r1 - r2) * glm::min(theta_step, t2 + theta_step - t1);
                    a2 = glm::min(radius_step, r1 - r2) * glm::min(theta_step, t1 - t2);
                    a3 = glm::min(radius_step, r2 + radius_step - r1) * glm::min(theta_step, t1 - t2);
                    n1 = 0;
                    n2 = 7;
                    n3 = 6;
                }

            }

            // area 1
            int i2 = n[n1];
            double dm = newSegment[i1].density * a1;
            double m1 = newSegment[i1].density * newSegment[i1].area;
            if (i2 < 0) {
                if (n1 == 0) {
                    *disk->getStellarMass() += dm;
                } else {
                    *disk->getEscapeMass() += dm;
                }
            } else {
                newSegment[i2].pr.push_back({newSegment[i1].vr, dm});
            }
            m1 = glm::max(0.0, m1 - dm);

            // area 2
            i2 = n[n2];
            dm = newSegment[i1].density * a2;
            if (i2 < 0) {
                if (n2 == 1 || n2 == 7) {
                    *disk->getStellarMass() += dm;
                } else {
                    *disk->getEscapeMass() += dm;
                }
            } else {

                newSegment[i2].pr.push_back({newSegment[i1].vr, dm / 2});
                newSegment[i2].pt.push_back({newSegment[i1].vt, dm / 2});

            }
            m1 = glm::max(0.0, m1 - dm);


            // area 3
            i2 = n[n3];
            dm = newSegment[i1].density * a3;
            newSegment[i2].pt.push_back({newSegment[i1].vt, dm});
            m1 = glm::max(0.0, m1 - dm);

            newSegment[i1].pr.push_back({newSegment[i1].vr, m1 / 2});
            newSegment[i1].pt.push_back({newSegment[i1].vt, m1 / 2});

        }
    }
    for (int i = 0; i < newSegment.size(); i++) {
        double p = 0.0, m = 0.0;
        for (auto &it : newSegment[i].pt) {
            p += it.m * it.v;
            m += it.m;
        }
        newSegment[i].vt = 0;
        if (m > 0) {
            newSegment[i].vt = p / m;
        }
        newSegment[i].density = m / newSegment[i].area;
        p = 0.0, m = 0.0;
        for (auto &it : newSegment[i].pr) {
            p += it.m * it.v;
            m += it.m;
        }
        newSegment[i].vr = 0;
        if (m > 0) {
            newSegment[i].vr = p / m;
        }
        newSegment[i].density += m / newSegment[i].area;
    }

}

void ApplyBruteForceGravity( Disk * disk ) {

    initialiseAcceleration(
            disk->getNewSegment(),
            disk->getSegment()
    );

    GravityProvider *bfgp = new PolarBruteForceSelfGravityProvider(
            disk->getNewSegment(),
            disk->getSegment(),
            disk->get_num_radial_cells(),
            disk->get_num_azimuthal_cells(),
            disk->getStellarMass()
    );

    bfgp->calculate();

    GravityProvider *sgp = new StellarGravityProvider(
            disk->getNewSegment(),
            disk->getSegment(),
            disk->get_num_radial_cells(),
            disk->get_num_azimuthal_cells(),
            disk->getStellarMass()
    );

    sgp->calculate();

    ApplyAcceleration(
            disk->getNewSegment()
    );

    MoveMass( disk );

}


/*
#define AFAC 0.25f
#define BFAC 0.5f

void ApplyXNearestNeighbourGravity() {

    std::vector<double> dd(num_points * num_points, 0.0f);
    int span = 3;

    for (std::size_t t = 0; t < num_points; t++) {

        for (std::size_t r = 0; r < num_points; r++) {

            std::size_t i1 = (t * num_points) + r;

            double dr = -1/(segment[i1].r*segment[i1].r);
            double dt = 0;

            for ( int st = -span; st < span; st++) {
                for ( int sr = -span; sr < span; sr++) {

                    std::size_t tr = max(0, glm::min((int) num_points, (int)r + sr));
                    std::size_t tt = ((t + num_points + st) % num_points);
                    std::size_t i2 = (tt * num_points) + tr;

                    if ( i2 != i1 ) {

                        double x = segment[i2].x - segment[i1].x;
                        double y = segment[i2].y - segment[i1].y;
                        double d = sqrt((x*x)*(y*y));
                        double f = d * d;
                        double nr = 0;
                        if (f > 0) {
                            nr = segment[i2].m / f;
                        }
                        double nt = 0;
                        if ( d > 0) {
                            nt = asin(y / d);
                        }
                        dt += nr * sin(nt);
                        dr += nr * cos(nt);
                    }
                }
            }

            int n[8] = { -1, -1, -1, -1, -1, -1, -1, -1};
            if ( r > 0 ) {
                n[0] = (t * num_points) + r - 1;
                n[1] = (((t+1) % num_points) * num_points) + r - 1;
                n[7] = (((t + num_points - 1) % num_points) * num_points) + r - 1;
            }
            if ( r < num_points - 1 ) {
                n[3] = (((t+1) % num_points) * num_points) + r + 1;
                n[4] = (t * num_points) + r + 1;
                n[5] = (((t + num_points - 1) % num_points) * num_points) + r + 1;

            }
            n[2] = (((t+1) % num_points) * num_points) + r;
            n[6] = (((t + num_points - 1) % num_points) * num_points) + r;


            double d = sqrt(dr*dr+dt*dt);
            d = min( segment[i1].m, d );
            dd[i1] -= d;

            double pdr;
            double pdt;
            double pc = (glm::abs(dr) + glm::abs(dt));
            if ( pc > 0 ) {
                pdr = d * glm::abs(dr) / (glm::abs(dr) + glm::abs(dt));
                pdt = d * glm::abs(dt) / (glm::abs(dr) + glm::abs(dt));
            } else {
                pdr = 0;
                pdt = 0;
            }

            if ( dr > 0 ) {
                if( n[3] >= 0 ) {
                    dd[n[3]] += AFAC * pdr;
                    dd[n[4]] += BFAC * pdr;
                    dd[n[5]] += AFAC * pdr;
                }
            } else {
                if( n[1] >= 0 ) {
                    dd[n[1]] += AFAC * pdr;
                    dd[n[0]] += BFAC * pdr;
                    dd[n[7]] += AFAC * pdr;
                }
            }

            if ( dt > 0 ) {
                dd[n[2]] += BFAC * pdt;
                if( n[1] >= 0 ) {

                    dd[n[1]] += AFAC * pdt;
                }
                if( n[3] >= 0 ) {
                    dd[n[3]] += AFAC * pdt;
                }
            } else {
                if( n[7] >= 0 ) {
                    dd[n[7]] += AFAC * pdt;
                }
                    if( n[5] >= 0 ) {
                        dd[n[5]] += AFAC * pdt;
                    }
                dd[n[6]] += BFAC * pdt;
            }

        }
    }

    for ( int i = 0; i < newSegment.size(); i++) {
        newSegment[i].m = max( 0.0, newSegment[i].m + 1e16*dd[i]);
    }
}

void ApplyNearestNeighbourGravity() {

    std::vector<double> dd(num_points * num_points, 0.0);


    for (std::size_t t = 0; t < num_points; t++) {

        for (std::size_t r = 1; r < num_points - 1; r++) {

            std::size_t i = (t * num_points) + r;
            std::size_t n[8];
            n[0] = (t * num_points) + r - 1;
            n[1] = (((t+1) % num_points) * num_points) + r - 1;
            n[2] = (((t+1) % num_points) * num_points) + r;
            n[3] = (((t+1) % num_points) * num_points) + r + 1;
            n[4] = (t * num_points) + r + 1;
            n[5] = (((t + num_points - 1) % num_points) * num_points) + r + 1;
            n[6] = (((t + num_points - 1) % num_points) * num_points) + r;
            n[7] = (((t + num_points - 1) % num_points) * num_points) + r - 1;

            // radial component
            double dr =   ( r2 * segment[n[3]].m + segment[n[4]].m + r2 * segment[n[5]].m )
                       - ( r2 * segment[n[1]].m + segment[n[0]].m + r2 * segment[n[7]].m )
                       - 10/(segment[i].r*segment[i].r);

            // angular component
            double dt =   ( r2 * segment[n[1]].m + segment[n[2]].m + r2 * segment[n[3]].m )
                       - ( r2 * segment[n[5]].m + segment[n[6]].m + r2 * segment[n[7]].m );

            double d = sqrt(dr*dr+dt*dt);
            d = min( segment[i].m, d );
            dd[i] -= d;

            double pdr;
            double pdt;
            double pc = (glm::abs(dr) + glm::abs(dt));
            if ( pc > 0 ) {
                pdr = d * glm::abs(dr) / (glm::abs(dr) + glm::abs(dt));
                pdt = d * glm::abs(dt) / (glm::abs(dr) + glm::abs(dt));
            } else {
                pdr = 0;
                pdt = 0;
            }
            if ( dr > 0 ) {
                dd[n[3]] += 0.25 * pdr;
                dd[n[4]] += 0.50 * pdr;
                dd[n[5]] += 0.25 * pdr;
            } else {
                dd[n[1]] += 0.25 * pdr;
                dd[n[0]] += 0.50 * pdr;
                dd[n[7]] += 0.25 * pdr;
            }

            if ( dt > 0 ) {
                dd[n[1]] += 0.25 * pdt;
                dd[n[2]] += 0.50 * pdt;
                dd[n[3]] += 0.25 * pdt;
            } else {
                dd[n[7]] += 0.25 * pdt;
                dd[n[6]] += 0.50 * pdt;
                dd[n[5]] += 0.25 * pdt;
            }

        }
    }

    for ( int i = 0; i < newSegment.size(); i++) {
        newSegment[i].m = max( 0.0, newSegment[i].m + dd[i]);
    }
}
*/