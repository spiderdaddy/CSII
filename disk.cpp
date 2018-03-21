//
// Created by lf on 20/03/18.
//

#include <cstdlib>
#include <cstring>
#include <cmath>
#include <iostream>
#include <vector>
#include <fstream>

#include <glm/common.hpp>
#include <algorithm>

#include "disk.h"
#include "gravity.h"


double stellar_mass = M_SUN;
double escape_mass = 0;

std::vector<SegmentVertices> segmentVertices;
std::vector<SegmentColours> segmentColours;
std::vector<Segment> segment;
std::vector<Segment> newSegment;

std::vector<SegmentVertices> getSegmentVertices() { return segmentVertices; }

std::vector<SegmentColours> getSegmentColours() { return segmentColours; }

std::vector<Segment> *getSegment() { return &segment; }

std::vector<Segment> *getNewSegment() { return &newSegment; }

double *getStellarMass() { return &stellar_mass; }

double *getEscapeMass() { return &escape_mass; }


void InitializeDisk() {

    std::vector<double> densities = getDensities();

    std::vector<double>::iterator dIterator = densities.begin();

    segmentVertices.reserve(NUM_RADIAL_CELLS * NUM_AZIMUTHAL_CELLS);
    segmentColours.reserve(NUM_RADIAL_CELLS * NUM_AZIMUTHAL_CELLS);
    segment.reserve(NUM_RADIAL_CELLS * NUM_AZIMUTHAL_CELLS);
    newSegment.reserve(NUM_RADIAL_CELLS * NUM_AZIMUTHAL_CELLS);

    srand(time(NULL));

    // populate points
    double theta_step = 2 * M_PI / NUM_AZIMUTHAL_CELLS;
    double theta_step_2 = theta_step / 2;

    double r_ratio = std::pow(OUTER_RADIUS/INNER_RADIUS, 1/(double)NUM_RADIAL_CELLS);
    double radius = INNER_RADIUS * sqrt(r_ratio);

    for (size_t r = 0; r < NUM_RADIAL_CELLS; r++) {

        double radius_lower = radius/sqrt(r_ratio);
        double radius_upper = radius*sqrt(r_ratio);

        for (size_t t = 0; t < NUM_AZIMUTHAL_CELLS; t++) {

            double theta = t * theta_step + theta_step_2;

            Segment s;
            s.pr.reserve(9);
            s.pt.reserve(9);
            s.r = radius;
            s.theta = theta;
            s.x = radius * cos(theta);
            s.y = radius * sin(theta);
            s.area = ( radius_upper - radius_lower ) * theta_step;
            s.m = *dIterator * s.area;
            dIterator++;
            s.vr = 0;
            s.vt = sqrt(G * stellar_mass / s.r);

            if (rand() % 2 == 1) {
                s.vt *= -1;
            }

            if (r > 0) {
                s.n[0] = ((r - 1) * NUM_AZIMUTHAL_CELLS) + t;
                s.n[1] = ((r - 1) * NUM_AZIMUTHAL_CELLS) + ((t + 1) % NUM_AZIMUTHAL_CELLS);
                s.n[7] = ((r - 1) * NUM_AZIMUTHAL_CELLS) + ((t + NUM_AZIMUTHAL_CELLS - 1) % NUM_AZIMUTHAL_CELLS);
            }
            if (r < NUM_RADIAL_CELLS - 1) {
                s.n[4] = ((r + 1) * NUM_AZIMUTHAL_CELLS) + t;
                s.n[3] = ((r + 1) * NUM_AZIMUTHAL_CELLS) + ((t + 1) % NUM_AZIMUTHAL_CELLS);
                s.n[5] = ((r + 1) * NUM_AZIMUTHAL_CELLS) + ((t + NUM_AZIMUTHAL_CELLS - 1) % NUM_AZIMUTHAL_CELLS);

            }
            s.n[2] = (r * NUM_AZIMUTHAL_CELLS) + ((t + 1) % NUM_AZIMUTHAL_CELLS);
            s.n[6] = (r * NUM_AZIMUTHAL_CELLS) + ((t + NUM_AZIMUTHAL_CELLS - 1) % NUM_AZIMUTHAL_CELLS);


            segment.push_back(s);
            newSegment.push_back(s);

            SegmentVertices sv;
            sv.v1.x = (float)(radius_lower * cos(theta - theta_step_2));
            sv.v1.y = (float)(radius_lower * sin(theta - theta_step_2));
            sv.v1.z = 0.0f;
            sv.v1.w = 1.0f;
            sv.v2.x = (float)(radius_upper * cos(theta - theta_step_2));
            sv.v2.y = (float)(radius_upper * sin(theta - theta_step_2));
            sv.v2.z = 0.0f;
            sv.v2.w = 1.0f;
            sv.v3.x = (float)(radius_upper * cos(theta + theta_step_2));
            sv.v3.y = (float)(radius_upper * sin(theta + theta_step_2));
            sv.v3.z = 0.0f;
            sv.v3.w = 1.0f;
            sv.v4.x = (float)(radius_lower  * cos(theta - theta_step_2));
            sv.v4.y = (float)(radius_lower  * sin(theta - theta_step_2));
            sv.v4.z = 0.0f;
            sv.v4.w = 1.0f;
            sv.v5.x = (float)(radius_lower  * cos(theta + theta_step_2));
            sv.v5.y = (float)(radius_lower  * sin(theta + theta_step_2));
            sv.v5.z = 0.0f;
            sv.v5.w = 1.0f;
            sv.v6.x = (float)(radius_upper * cos(theta + theta_step_2));
            sv.v6.y = (float)(radius_upper * sin(theta + theta_step_2));
            sv.v6.z = 0.0f;
            sv.v6.w = 1.0f;
            segmentVertices.push_back(sv);

            SegmentColours sc;
            sc.c1.r = s.m;
            sc.c1.g = sc.c1.r;
            sc.c1.b = sc.c1.r;
            sc.c1.a = 1.0f;
            sc.c2 = sc.c1;
            sc.c3 = sc.c1;
            sc.c4 = sc.c1;
            sc.c5 = sc.c1;
            sc.c6 = sc.c1;
            segmentColours.push_back(sc);
        }

        radius *= r_ratio;
    }

    fprintf(
            stderr,
            "INFO: segmentVertices.size(): %d\n",
            (unsigned) segmentVertices.size()
    );

    //PrintPoints();
    CalcSystemMass();
}


void PrintPoints() {
    for (int i = 0; i < segmentVertices.size(); i++) {
        SegmentVertices pt = segmentVertices[i];
        fprintf(stdout,
                "%d:([%f,%f],[%f,%f],[%f,%f]) ",
                i, pt.v1.x, pt.v1.y, pt.v2.x, pt.v2.y, pt.v3.x, pt.v3.y);
    }
    fprintf(stdout, "\n");
}

void CalcSystemMass() {
    double disk_mass = 0;
    double pr = 0;
    double pt = 0;
    for (std::size_t i = 0; i < segment.size(); i++) {
        disk_mass += segment[i].m;
        pr += segment[i].m * segment[i].vr;
        pt += segment[i].m * segment[i].vt;
    }
    fprintf(stdout,
            "INFO: Stellar Mass: %.8e Disk Mass: %.8e : Escape Mass: %.8e Total mass: %.8e, pr = %.8e, pt = %.8e\n",
            stellar_mass, disk_mass, escape_mass, stellar_mass + disk_mass, pr, pt);
}

void MapSegmentToColor() {

    double d_max = 250;

    for (int i = 0; i < segmentColours.size(); i++) {
        SegmentColours *scp = &segmentColours[i];
        scp->c1.r = std::min( std::abs( std::min((float)(segment[i].m / segment[i].area ), 0.0f) ), 1.0f );
        scp->c1.g = std::min( std::abs( std::max((float)(segment[i].m / segment[i].area / 3.0f), 0.0f) ), 1.0f );
        scp->c1.b = 0.15f;
        scp->c2 = scp->c1;
        scp->c3 = scp->c1;
        scp->c4 = scp->c1;
        scp->c5 = scp->c1;
        scp->c6 = scp->c1;
    }
}

void swapSegments() {

    segment = newSegment;
}


using namespace std;

double minimum;
double maximum;
double span;
std::vector<double> getDensities() {

    std::vector<double> densities;
    string line;
    ifstream myfile("/data/UZH/CSII/data1/density.data");
    if (myfile.is_open()) {

        getline( myfile, line );
        while (myfile.good()) {
            if( line.size() > 0 ) {
                double density = stod(line);
                densities.push_back(density);
            }
            getline( myfile, line );
        }
        myfile.close();
    } else cout << "Unable to open file";

    minimum = *std::min_element(densities.begin(), densities.end());
    maximum = *std::max_element(densities.begin(), densities.end());
    fprintf(stdout,
            "INFO: Density: Minimum %.8e Maximum: %.8e\n", minimum, maximum);


    /*
    span = maximum - minimum;
    minimum = abs(minimum);
    std::for_each(densities.begin(), densities.end(), [](double &el){el = (el+minimum)/span ; });
    */
    return densities;
}