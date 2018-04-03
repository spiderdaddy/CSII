//
// Created by lf on 21/03/18.
//

#include <chrono>
#include "../disk.h"

#include "GravityProvider.h"


GravityProvider::GravityProvider( Disk * d ) {
    disk = d;
    pNewSegment = disk->getNewSegment();
    pSegment = disk->getSegment();
    num_radial_cells = disk->get_num_radial_cells();
    num_azimuthal_cells = disk->get_num_azimuthal_cells();
    stellar_mass = disk->getStellarMass();
}

void GravityProvider::calculate() {

}

void GravityProvider::startTimer() {
    startTime = std::chrono::system_clock::now();
}

void GravityProvider::stopTimer() {
    endTime = std::chrono::system_clock::now();

}

long GravityProvider::getTime() {
    return std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();
}
