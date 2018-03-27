//
// Created by lf on 20/03/18.
//
#ifndef CSII_DISK_H
#define CSII_DISK_H

#include <string>
#include <vector>

#define AU 1.49e08
#define INNER_RADIUS 0.5
#define OUTER_RADIUS 1.5

#define G 6.6e-11
#define M_EARTH   1e22
#define M_SUN     2e30
#define M_JUPITER 1e27
#define STAR
#define SELF

class Disk {

public:

    Disk(unsigned int num_raidal_cells, unsigned int num_azimuthal_cells, std::string filename );

    struct XYZW_GL {
        float x, y, z, w;
    };

    struct RGBA_GL {
        float r, g, b, a;
    };

    struct SegmentVertices {
        XYZW_GL v1;
        XYZW_GL v2;
        XYZW_GL v3;
        XYZW_GL v4;
        XYZW_GL v5;
        XYZW_GL v6;

    };
    struct SegmentColours {
        RGBA_GL c1;
        RGBA_GL c2;
        RGBA_GL c3;
        RGBA_GL c4;
        RGBA_GL c5;
        RGBA_GL c6;

    };

    struct Polar {
        double r;
        double t;
    };

    struct P {
        double v, m;
    };

    struct Segment {
        double r, theta, x, y;
        double vr, vt;
        double ar, at;
        double area;
        double density;
        long n[8] = {-1, -1, -1, -1, -1, -1, -1, -1}; //nearest neighbours
        std::vector<P> pr;
        std::vector<P> pt;
    };

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

    std::vector<double> loadDensities(std::string filename);

    unsigned int get_num_radial_cells();

    unsigned int get_num_azimuthal_cells();

    double get_stellar_mass();

    double get_escape_mass();

private:
    unsigned int num_radial_cells;
    unsigned int num_azimuthal_cells;
    double stellar_mass = M_SUN;
    double escape_mass = 0;


};

#endif //CSII_DISK_H
