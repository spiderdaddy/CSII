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


using namespace std;


std::vector<Disk::SegmentVertices> Disk::getSegmentVertices() { return segmentVertices; }

std::vector<Disk::SegmentVertices> Disk::getGravityVertices() { return gravityVertices; }

std::vector<Disk::SegmentColours> Disk::getSegmentColours() { return segmentColours; }

std::vector<Disk::SegmentColours> Disk::getGravityColours() { return gravityColours; }

std::vector<Disk::Segment> *Disk::getSegment() { return &segment; }

std::vector<Disk::Segment> *Disk::getNewSegment() { return &newSegment; }

double *Disk::getStellarMass() { return &stellar_mass; }

double *Disk::getEscapeMass() { return &escape_mass; }

int Disk::getCellIndex(int r, int theta) {
    return ( r * num_azimuthal_cells ) + theta;
}

double Disk::getDensity( int r, int theta ) {
    return segment[getCellIndex(r, theta)].density;
}

double Disk::getArea( int r, int theta ) {
    return segment[getCellIndex(r, theta)].area;
}

double Disk::getRLower( int r, int theta ) {
    return segment[getCellIndex(r, theta)].r / sqrt_r_ratio;
}

double Disk::getRUpper( int r, int theta ) {
    return segment[getCellIndex(r, theta)].r * sqrt_r_ratio;
}

double Disk::getThetaLower( int r, int theta ) {
    return segment[getCellIndex(r, theta)].theta - theta_step_2;
}

double Disk::getThetaUpper( int r, int theta ) {
    return segment[getCellIndex(r, theta)].theta + theta_step_2;
}

void Disk::setTreeNode( int r, int theta, QTNode * node) {
    segment[getCellIndex(r, theta)].node = node;
}


Disk::Disk( unsigned int num_r, unsigned int num_theta, std::string filename ) {

    num_radial_cells = num_r;
    num_azimuthal_cells = num_theta;

    std::vector<double> densities = loadDensities(filename);

    std::vector<double>::iterator dIterator = densities.begin();

    segmentVertices.reserve(num_radial_cells * num_azimuthal_cells);
    segmentColours.reserve(num_radial_cells * num_azimuthal_cells);
    gravityVertices.reserve(num_radial_cells * num_azimuthal_cells);
    gravityColours.reserve(num_radial_cells * num_azimuthal_cells);
    segment.reserve(num_radial_cells * num_azimuthal_cells);
    newSegment.reserve(num_radial_cells * num_azimuthal_cells);

    srand(time(NULL));

    // populate points
    theta_step = 2 * M_PI / num_azimuthal_cells;
    theta_step_2 = theta_step / 2;

    r_ratio = std::pow(OUTER_RADIUS/INNER_RADIUS, 1/(double)num_radial_cells);
    sqrt_r_ratio = sqrt(r_ratio);
    double radius = INNER_RADIUS * sqrt(r_ratio);

    for (size_t r = 0; r < num_radial_cells; r++) {

        double radius_lower = radius/sqrt_r_ratio;
        double radius_upper = radius*sqrt_r_ratio;

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
            s.vt = sqrt(1 / s.r) * theta_step;

            if (r > 0) {
                s.neighbour[0] = ((r - 1) * num_azimuthal_cells) + t;
                s.neighbour[1] = ((r - 1) * num_azimuthal_cells) + ((t + 1) % num_azimuthal_cells);
                s.neighbour[7] = ((r - 1) * num_azimuthal_cells) + ((t + num_azimuthal_cells - 1) % num_azimuthal_cells);
            }
            if (r < num_radial_cells - 1) {
                s.neighbour[4] = ((r + 1) * num_azimuthal_cells) + t;
                s.neighbour[3] = ((r + 1) * num_azimuthal_cells) + ((t + 1) % num_azimuthal_cells);
                s.neighbour[5] = ((r + 1) * num_azimuthal_cells) + ((t + num_azimuthal_cells - 1) % num_azimuthal_cells);

            }
            s.neighbour[2] = (r * num_azimuthal_cells) + ((t + 1) % num_azimuthal_cells);
            s.neighbour[6] = (r * num_azimuthal_cells) + ((t + num_azimuthal_cells - 1) % num_azimuthal_cells);


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
            sv.v1.x -= 1.1 * OUTER_RADIUS;
            sv.v2.x -= 1.1 * OUTER_RADIUS;
            sv.v3.x -= 1.1 * OUTER_RADIUS;
            sv.v4.x -= 1.1 * OUTER_RADIUS;
            sv.v5.x -= 1.1 * OUTER_RADIUS;
            sv.v6.x -= 1.1 * OUTER_RADIUS;

            segmentVertices.push_back(sv);
            sv.v1.x += 2.2 * OUTER_RADIUS;
            sv.v2.x += 2.2 * OUTER_RADIUS;
            sv.v3.x += 2.2 * OUTER_RADIUS;
            sv.v4.x += 2.2 * OUTER_RADIUS;
            sv.v5.x += 2.2 * OUTER_RADIUS;
            sv.v6.x += 2.2 * OUTER_RADIUS;
            gravityVertices.push_back(sv);

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
            gravityColours.push_back(sc);
        }

        radius *= r_ratio;
    }

    qt = new QuadTree(num_radial_cells, num_azimuthal_cells);
    qt->printNodes();

    fprintf(
            stdout,
            "INFO: segmentVertices.size(): %d\n",
            (unsigned) segmentVertices.size()
    );

    MapSegmentToColor();
    MapGravityToColor();
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

    double d_max_pos = 0;
    double d_max_neg = 0;
    for (int i = 0; i < segmentColours.size(); i++) {
        d_max_pos = std::max(d_max_pos, segment[i].density);
        d_max_neg = std::min(d_max_neg, segment[i].density);
    }

    for (int i = 0; i < segmentColours.size(); i++) {
        SegmentColours *scp = &segmentColours[i];
        scp->c1.r = std::min( std::abs( std::min((float)(segment[i].density / abs(d_max_neg)), 0.0f) ), 1.0f );
        scp->c1.g = std::min( std::abs( std::max((float)(segment[i].density / abs(d_max_pos)), 0.0f) ), 1.0f );
        scp->c1.b = 0.15f;
        scp->c2 = scp->c1;
        scp->c3 = scp->c1;
        scp->c4 = scp->c1;
        scp->c5 = scp->c1;
        scp->c6 = scp->c1;
    }
}

/*
 * https://stackoverflow.com/questions/37876316/map-value-range-to-rainbow-colormap
 * Please try this function, which has linear interpolation and wraps around to red->green->blue->red.
 * np is your maximum value (a) and p is the input value (v).
 * You could get it to stop at violet by increasing np a bit so that p is always less than np.
 */

void getcolor(double p, double np, float&r, float&g, float&b) {
    double inc = 6.0 / np;
    double x = p * inc;
    r = 0.0f; g = 0.0f; b = 0.0f;
    if ((0 <= x && x <= 1) || (5 <= x && x <= 6)) r = 1.0f;
    else if (4 <= x && x <= 5) r = x - 4;
    else if (1 <= x && x <= 2) r = 1.0f - (x - 1);
    if (1 <= x && x <= 3) g = 1.0f;
    else if (0 <= x && x <= 1) g = x - 0;
    else if (3 <= x && x <= 4) g = 1.0f - (x - 3);
    if (3 <= x && x <= 5) b = 1.0f;
    else if (2 <= x && x <= 3) b = x - 2;
    else if (5 <= x && x <= 6) b = 1.0f - (x - 5);
}

void Disk::MapGravityToColor() {

    double d_max = 0;

    for (int i = 0; i < gravityColours.size(); i++) {
        d_max = max(d_max, sqrt(segment[i].at*segment[i].at + segment[i].ar*segment[i].ar));
    }

    for (int i = 0; i < gravityColours.size(); i++) {
        SegmentColours *scp = &gravityColours[i];
        double H = sqrt(segment[i].at*segment[i].at + segment[i].ar*segment[i].ar);
        double theta = 0;
        if ( segment[i].ar > 0 ) theta = atan(segment[i].at/segment[i].ar);
        else theta = atan(segment[i].at/segment[i].ar) + M_PI;
//        fprintf(
//                stdout,
//                "INFO: %f %f %0.18f %0.18f\r",
//                segment[i].r, segment[i].theta, H, theta
//        );
        theta += segment[i].theta;
        if (theta < 0 ) theta = (2.0 * M_PI) + theta;
        if (theta > 2.0 * M_PI) theta = theta - (2.0 * M_PI);
        getcolor( theta, (2.0 * M_PI), scp->c1.r, scp->c1.g, scp->c1.b );
        double N = H / d_max;
        scp->c1.r *= N;
        scp->c1.g *= N;
        scp->c1.b *= N;
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

std::vector<double> Disk::loadDensities( string filename ) {

    std::vector<double> densities;
    string line;
    ifstream myfile(filename);
    if (myfile.is_open()) {

        getline( myfile, line );
        while (myfile.good()) {
            if( line.size() > 0 ) {
                double density = stod(line);
                //density = 0.01;
                densities.push_back(density);
            }
            getline( myfile, line );
        }
        myfile.close();
    } else cout << "Unable to open file";

    minimum_density = *std::min_element(densities.begin(), densities.end());
    maximum_density = *std::max_element(densities.begin(), densities.end());
    fprintf(stdout,
            "INFO: Density: Minimum %.8e Maximum: %.8e\n", minimum_density, maximum_density);


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

QuadTree *Disk::getQuadTree() const {
    return qt;
}


QTNode *Disk::calcTreeValues(QTNode *parent, QTNode *node) {

    // If this node is a single cell, then just set the density and return
    if ( (node->r_start == node->r_end) &&
         (node->t_start == node->t_end) ) {
        node->density = segment[getCellIndex( node->r_start, node->t_start )].density;
        node->area = segment[getCellIndex( node->r_start, node->t_start )].area;
        node->r = segment[getCellIndex( node->r_start, node->t_start )].r;
        node->theta = segment[getCellIndex( node->r_start, node->t_start )].theta;
        segment[getCellIndex( node->r_start, node->t_start )].node =  node;
    } else {

        int cells = 0;
        node->area = 0;
        node->density = 0;

        for (int i = 0; i < 4; i++) {
            if(node->leaf[i] != nullptr) {
                cells ++;
                calcTreeValues(node, node->leaf[i]);
                node->area += node->leaf[i]->area;
                node->density += node->leaf[i]->density;
            }
        }

        node->density /= cells;

        node->r = ( segment[getCellIndex( node->r_start, node->t_start )].r + segment[getCellIndex( node->r_end, node->t_end )].r ) / 2.0;
        node->theta = ( segment[getCellIndex( node->r_start, node->t_start )].theta + segment[getCellIndex( node->r_end, node->t_end )].theta ) / 2.0;

    }
    return node;
}


