//
// Created by lf on 21/03/18.
//

#ifndef CSII_GRAVITYPROVIDER_H
#define CSII_GRAVITYPROVIDER_H

#include <chrono>
#include "../disk.h"


class GravityProvider {

public:

    GravityProvider()= default;

    GravityProvider( Disk * disk);

    virtual void calculate();
    long getTime();
    void setResolution(int);
    void setDepth(int);


protected:
    Disk *disk;
    std::vector<Disk::Segment> *pNewSegment;
    std::vector<Disk::Segment> *pSegment;
    int num_radial_cells;
    int num_azimuthal_cells;
    double * stellar_mass;
    void startTimer();
    void stopTimer();
    int resolution = 0;
    int depth = 0;

private:
    std::chrono::system_clock::time_point startTime;
    std::chrono::system_clock::time_point endTime;
};


#endif //CSII_GRAVITYPROVIDER_H
