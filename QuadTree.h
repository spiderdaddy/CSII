//
// Created by lf on 27/03/18.
//

#ifndef CSII_QUADTREE_H
#define CSII_QUADTREE_H

#include <iostream>

using namespace std;

struct QTNode {

public:
    QTNode(QTNode *p, int rs, int re, int ts, int te, unsigned tree_level );

    QTNode *parent;
    QTNode *leaf[4];
    int r_start;
    int r_end;
    int t_start;
    int t_end;
    unsigned tree_level;
    double area;
    double density;
    double r;
    double theta;
};

class QuadTree {
public:

    QuadTree( unsigned num_radial_points, unsigned num_azimuthal_points);

    void printNodes();

private:
    QTNode* createNode(QTNode *parent, int r_start, int r_end, int t_start, int t_end, unsigned tree_level );
    unsigned max_level = 0;
    unsigned num_radial_points;
    unsigned num_azimuthal_points;
    QTNode *head;
public:
    QTNode *getHead() const;

};


#endif //CSII_QUADTREE_H
