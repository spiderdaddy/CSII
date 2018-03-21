//
// Created by lf on 20/03/18.
//
#ifndef CSII_DISK_H
#define CSII_DISK_H

#include <vector>

#define NUM_RADIAL_CELLS   128 //128
#define NUM_AZIUTHAL_CELLS 256 //256


#define AU 1.49e08
#define INNER_RADIUS (0.39*AU)
#define OUTER_RADIUS (40.00*AU)

struct XY {
    float x, y, z, w;
};
struct RGBA {
    float r, g, b, a;
};
struct SegmentVertices {
    XY v1;
    XY v2;
    XY v3;
    XY v4;
    XY v5;
    XY v6;

};
struct SegmentColours {
    RGBA c1;
    RGBA c2;
    RGBA c3;
    RGBA c4;
    RGBA c5;
    RGBA c6;

};
struct P {
    double v, m;
};
struct Segment {
    double r, theta, x, y;
    double vr, vt;
    double ar, at;
    double area;
    double m;
    long n[8] = { -1, -1, -1, -1, -1, -1, -1, -1}; //nearest neighbours
    std::vector<P> pr;
    std::vector<P> pt;
};

void InitializeDisk();

void PrintPoints();

void CalcSystemMass();

void MapSegmentToColor();

std::vector<SegmentVertices> getSegmentVertices();
std::vector<SegmentColours> getSegmentColours();
std::vector<Segment> *getSegment();
std::vector<Segment> *getNewSegment();
double *getStellarMass();
double *getEscapeMass();
void swapSegments();

#endif //CSII_DISK_H
