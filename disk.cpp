//
// Created by lf on 20/03/18.
//

#include <cstring>
#include <cmath>
#include <iostream>

#include <vector>
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
std::vector<Segment> getSegment(){ return segment; }
std::vector<Segment> getNewSegment(){return newSegment; }
double *getStellarMass() { return &stellar_mass; }
double *getEscapeMass() { return &escape_mass;}
void swapSegments() { segment = newSegment; }


void InitializeDisk() {
    segmentVertices.reserve(NUM_RADIAL_CELLS * NUM_AZIUTHAL_CELLS);
    segmentColours.reserve(NUM_RADIAL_CELLS * NUM_AZIUTHAL_CELLS);
    segment.reserve(NUM_RADIAL_CELLS * NUM_AZIUTHAL_CELLS);
    newSegment.reserve(NUM_RADIAL_CELLS * NUM_AZIUTHAL_CELLS);

    srand (time(NULL));

    // populate points
    double theta_step = 2 * M_PI / NUM_AZIUTHAL_CELLS;
    double theta_step_2 = theta_step / 2;
    double radius_step = ( OUTER_RADIUS - INNER_RADIUS ) / NUM_RADIAL_CELLS;
    double radius_step_2 = radius_step / 2;

    for (size_t t = 0; t < NUM_AZIUTHAL_CELLS; t++) {

        double theta = t * theta_step + theta_step_2;

        for (size_t r = 0; r < NUM_RADIAL_CELLS; r++) {

            double radius = INNER_RADIUS + (r * radius_step) + radius_step_2;

            Segment s;
            s.pr.reserve(9);
            s.pt.reserve(9);
            s.r = radius;
            s.theta = theta;
            s.x = radius * cos(theta);
            s.y = radius * sin(theta);
            s.area = 2.0 * radius * radius_step * theta_step;
            s.m = (double)(rand()%50) * s.area;
            s.vr = 0;
            s.vt = sqrt(G * stellar_mass / s.r);

            if (rand()%2 == 1) {
                s.vt *= -1;
            }

            if ( r > 0 ) {
                s.n[0] = (t * NUM_RADIAL_CELLS) + r - 1;
                s.n[1] = (((t+1) % NUM_RADIAL_CELLS) * NUM_RADIAL_CELLS) + r - 1;
                s.n[7] = (((t + NUM_RADIAL_CELLS - 1) % NUM_RADIAL_CELLS) * NUM_RADIAL_CELLS) + r - 1;
            }
            if ( r < NUM_RADIAL_CELLS - 1 ) {
                s.n[3] = (((t+1) % NUM_RADIAL_CELLS) * NUM_RADIAL_CELLS) + r + 1;
                s.n[4] = (t * NUM_RADIAL_CELLS) + r + 1;
                s.n[5] = (((t + NUM_RADIAL_CELLS - 1) % NUM_RADIAL_CELLS) * NUM_RADIAL_CELLS) + r + 1;

            }
            s.n[2] = (((t+1) % NUM_RADIAL_CELLS) * NUM_RADIAL_CELLS) + r;
            s.n[6] = (((t + NUM_RADIAL_CELLS - 1) % NUM_RADIAL_CELLS) * NUM_RADIAL_CELLS) + r;


            segment.push_back(s);
            newSegment.push_back(s);

            SegmentVertices sv;
            sv.v1.x = (radius - radius_step_2) * cos(theta - theta_step_2);
            sv.v1.y = (radius - radius_step_2) * sin(theta - theta_step_2);
            sv.v1.z = 0.0f;
            sv.v1.w = 1.0f;
            sv.v2.x = (radius + radius_step_2) * cos(theta - theta_step_2);
            sv.v2.y = (radius + radius_step_2) * sin(theta - theta_step_2);
            sv.v2.z = 0.0f;
            sv.v2.w = 1.0f;
            sv.v3.x = (radius + radius_step_2) * cos(theta + theta_step_2);
            sv.v3.y = (radius + radius_step_2) * sin(theta + theta_step_2);
            sv.v3.z = 0.0f;
            sv.v3.w = 1.0f;
            sv.v4.x = (radius - radius_step_2) * cos(theta - theta_step_2);
            sv.v4.y = (radius - radius_step_2) * sin(theta - theta_step_2);
            sv.v4.z = 0.0f;
            sv.v4.w = 1.0f;
            sv.v5.x = (radius - radius_step_2) * cos(theta + theta_step_2);
            sv.v5.y = (radius - radius_step_2) * sin(theta + theta_step_2);
            sv.v5.z = 0.0f;
            sv.v5.w = 1.0f;
            sv.v6.x = (radius + radius_step_2) * cos(theta + theta_step_2);
            sv.v6.y = (radius + radius_step_2) * sin(theta + theta_step_2);
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
    }

    fprintf(
            stderr,
            "INFO: segmentVertices.size(): %d\n",
            (unsigned) segmentVertices.size()
    );

    PrintPoints();
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
    for (std::size_t i = 0; i < segment.size(); i ++ ) {
        disk_mass += segment[i].m;
        pr += segment[i].m * segment[i].vr;
        pt += segment[i].m * segment[i].vt;
    }
    fprintf(stdout, "INFO: Stellar Mass: %.8e Disk Mass: %.8e : Escape Mass: %.8e Total mass: %.8e, pr = %.8e, pt = %.8e\n",
            stellar_mass, disk_mass, escape_mass, stellar_mass+disk_mass, pr, pt );
}
