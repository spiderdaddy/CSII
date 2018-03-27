//
// Created by lf on 20/03/18.
//

#include <cstdlib>
#include <cstring>
#include <cmath>
#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>

#include <glm/common.hpp>


#include "disk.h"


std::vector<Disk::SegmentVertices> segmentVertices;
std::vector<Disk::SegmentColours> segmentColours;
std::vector<Disk::Segment> segment;
std::vector<Disk::Segment> newSegment;

std::vector<Disk::SegmentVertices> Disk::getSegmentVertices() { return segmentVertices; }

std::vector<Disk::SegmentColours> Disk::getSegmentColours() { return segmentColours; }

std::vector<Disk::Segment> *Disk::getSegment() { return &segment; }

std::vector<Disk::Segment> *Disk::getNewSegment() { return &newSegment; }

double *Disk::getStellarMass() { return &stellar_mass; }

double *Disk::getEscapeMass() { return &escape_mass; }


Disk::Disk( unsigned int num_r, unsigned int num_theta, std::string filename ) {

    num_radial_cells = num_r;
    num_azimuthal_cells = num_theta;

    std::vector<double> densities = loadDensities(filename);

    std::vector<double>::iterator dIterator = densities.begin();

    segmentVertices.reserve(num_radial_cells * num_azimuthal_cells);
    segmentColours.reserve(num_radial_cells * num_azimuthal_cells);
    segment.reserve(num_radial_cells * num_azimuthal_cells);
    newSegment.reserve(num_radial_cells * num_azimuthal_cells);

    srand(time(NULL));

    // populate points
    double theta_step = 2 * M_PI / num_azimuthal_cells;
    double theta_step_2 = theta_step / 2;

    double r_ratio = std::pow(OUTER_RADIUS/INNER_RADIUS, 1/(double)num_radial_cells);
    double radius = INNER_RADIUS * sqrt(r_ratio);

    for (size_t r = 0; r < num_radial_cells; r++) {

        double radius_lower = radius/sqrt(r_ratio);
        double radius_upper = radius*sqrt(r_ratio);

        for (size_t t = 0; t < num_azimuthal_cells; t++) {

            double theta = t * theta_step + theta_step_2;

            Segment s;
            s.pr.reserve(9);
            s.pt.reserve(9);
            s.r = radius;
            s.theta = theta;
            s.x = radius * cos(theta);
            s.y = radius * sin(theta);
            s.area = ( radius_upper - radius_lower ) * theta_step;
            s.density = *dIterator;
            dIterator++;
            s.vr = 0;
            s.vt = sqrt(G * stellar_mass / s.r);

            if (rand() % 2 == 1) {
                s.vt *= -1;
            }

            if (r > 0) {
                s.n[0] = ((r - 1) * num_azimuthal_cells) + t;
                s.n[1] = ((r - 1) * num_azimuthal_cells) + ((t + 1) % num_azimuthal_cells);
                s.n[7] = ((r - 1) * num_azimuthal_cells) + ((t + num_azimuthal_cells - 1) % num_azimuthal_cells);
            }
            if (r < num_radial_cells - 1) {
                s.n[4] = ((r + 1) * num_azimuthal_cells) + t;
                s.n[3] = ((r + 1) * num_azimuthal_cells) + ((t + 1) % num_azimuthal_cells);
                s.n[5] = ((r + 1) * num_azimuthal_cells) + ((t + num_azimuthal_cells - 1) % num_azimuthal_cells);

            }
            s.n[2] = (r * num_azimuthal_cells) + ((t + 1) % num_azimuthal_cells);
            s.n[6] = (r * num_azimuthal_cells) + ((t + num_azimuthal_cells - 1) % num_azimuthal_cells);


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
            sc.c1.r = 0;
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
            stdout,
            "INFO: segmentVertices.size(): %d\n",
            (unsigned) segmentVertices.size()
    );

    MapSegmentToColor();
    //PrintPoints();
    CalcSystemMass();
}


void Disk::PrintPoints() {
    for (int i = 0; i < segmentVertices.size(); i++) {
        SegmentVertices pt = segmentVertices[i];
        fprintf(stdout,
                "%d:([%f,%f],[%f,%f],[%f,%f]) ",
                i, pt.v1.x, pt.v1.y, pt.v2.x, pt.v2.y, pt.v3.x, pt.v3.y);
    }
    fprintf(stdout, "\n");
}

void Disk::CalcSystemMass() {
    double disk_mass = 0;
    double pr = 0;
    double pt = 0;
    for (std::size_t i = 0; i < segment.size(); i++) {
        double m = segment[i].density * segment[i].area;
        disk_mass += m;
        pr += m * segment[i].vr;
        pt += m * segment[i].vt;
    }
    fprintf(stdout,
            "INFO: Stellar Mass: %.8e Disk Mass: %.8e : Escape Mass: %.8e Total mass: %.8e, pr = %.8e, pt = %.8e\n",
            stellar_mass, disk_mass, escape_mass, stellar_mass + disk_mass, pr, pt);
}

void Disk::MapSegmentToColor() {

    double d_max = 250;

    for (int i = 0; i < segmentColours.size(); i++) {
        SegmentColours *scp = &segmentColours[i];
        scp->c1.r = std::min( std::abs( std::min((float)segment[i].density, 0.0f) ), 1.0f );
        scp->c1.g = std::min( std::abs( std::max((float)segment[i].density / 3.0f, 0.0f) ), 1.0f );
        scp->c1.b = 0.15f;
        scp->c2 = scp->c1;
        scp->c3 = scp->c1;
        scp->c4 = scp->c1;
        scp->c5 = scp->c1;
        scp->c6 = scp->c1;
    }
}

void Disk::swapSegments() {

    segment = newSegment;
}


using namespace std;

double minimum;
double maximum;
double span;
std::vector<double> Disk::loadDensities( string filename ) {

    std::vector<double> densities;
    string line;
    ifstream myfile(filename);
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

unsigned int Disk::get_num_radial_cells() {
    return num_radial_cells;
}

unsigned int Disk::get_num_azimuthal_cells() {
    return num_azimuthal_cells;
}

double Disk::get_stellar_mass() {
    return stellar_mass;
}

double Disk::get_escape_mass() {
    return escape_mass;
}


