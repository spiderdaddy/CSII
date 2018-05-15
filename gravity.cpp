//
// Created by lf on 20/03/18.
//

#include <vector>
#include <cstdlib>
#include <algorithm>
#include <fstream>

#include <glm/glm.hpp>

#include "disk.h"
#include "gravity.h"
#include "CartesianBruteForceSelfGravityProvider.h"
#include "PolarBruteForceSelfGravityProvider.h"
#include "SimplePolarTreeSelfGravityProvider.h"
#include "ExclusionPolarTreeSelfGravityProvider.h"
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

//    std::for_each(segment.begin(), segment.end(), [](Disk::Segment &el) {
//        el.ar = 0;
//        el.at = 0;
//    });
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

            long *n = &newSegment[i1].neighbour[0];

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

void ApplyGravity(Disk *disk, int resolution, string data_name) {

    initialiseAcceleration(
            disk->getNewSegment(),
            disk->getSegment()
    );

    GravityProvider *gp = new ExclusionPolarTreeSelfGravityProvider( disk );
//    GravityProvider *gp = new PolarBruteForceSelfGravityProvider( disk );

    gp->setResolution(resolution);
    gp->calculate();

    string filename(100, 0);
    sprintf( &filename[0], "result/PolarTree-res-%s-%d.csv", &data_name[0], resolution);
    disk->saveGravities( filename );

    sprintf( &filename[0], "result/PolarTree-times-%s.csv", &data_name[0]);
    ofstream myfile;
    myfile.open(filename, std::ofstream::out | std::ofstream::app );
    if (myfile.is_open()) {
        myfile << resolution << "," << gp->getTime() << "\n";
        myfile.close();
    } else cout << "Unable to open file";


//    GravityProvider *gp = new PolarBruteForceSelfGravityProvider( disk );
//
//    gp->calculate();

    fprintf(
            stdout,
            "INFO: time required(s)=%f\n",
            gp->getTime() / 1000.0
    );

/*
    GravityProvider *sgp = new StellarGravityProvider( disk );

    sgp->calculate();
*/
    /*
    ApplyAcceleration(
            disk->getNewSegment()
    );

    MoveMass( disk );
    */
}


