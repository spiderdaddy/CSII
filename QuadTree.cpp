//
// Created by lf on 27/03/18.
//

#include <cstdlib>
#include <cstring>
#include <cmath>
#include <iostream>
#include <vector>
#include <bits/stdc++.h>

#include "QuadTree.h"

using namespace std;


QTNode::QTNode(QTNode *p, int rs, int re, int ts, int te, unsigned tree_l) {
    parent = p;
    leaf[0] = (QTNode *) nullptr;
    leaf[1] = (QTNode *) nullptr;
    leaf[2] = (QTNode *) nullptr;
    leaf[3] = (QTNode *) nullptr;
    r_start = rs;
    r_end = re;
    t_start = ts;
    t_end = te;
    tree_level = tree_l;
    area = density = r = theta = 0;
}

ostream &operator<<(ostream &os, const QTNode *node) {
    os << (void *) node << ": "
       << node->tree_level;

    QTNode* node_counter = (QTNode*)node;
    while( node_counter->parent != nullptr ) {
        os << "|   ";
        node_counter = node_counter->parent;
    }

    os << (void *) node->parent << ", ("
       << (void *) node->leaf[0] << ", " << (void *) node->leaf[1] << ", " << (void *) node->leaf[2] << ", "
       << (void *) node->leaf[3] << "), ("
       << node->r_start << ", " << node->r_end << "), ("
       << node->t_start << ", " << node->t_end << "), ("
       << node->area << ", " << node->density << ", " << node->r << ", " << node->theta << "), (";

    for (QTNode *neighbour : node->neighbour) {
        cout << (void *) neighbour << ",";
    }
    cout << ")";
    return os;
}

int QuadTree::calcIndex( int level, int r, int theta ) {

    r = r >> level;
    theta = theta >> level;
    int n = num_azimuthal_points >> level;
    return r*n+theta;
}

QuadTree::QuadTree(unsigned num_r, unsigned num_a) {
    num_radial_points = num_r;
    num_azimuthal_points = num_a;

    // Find out the maximum number of levels
    unsigned largest_axis = max(num_radial_points, num_azimuthal_points);
    while (pow(2, max_level) < largest_axis) {
        max_level++;
    }

    cout << "INFO: QuadTree: r=" << num_radial_points << " theta=" << num_azimuthal_points << " max_level=" << max_level << "\n";

    // Create the list of vectors which will contain the nodes for each level
    levelVector = std::vector<std::vector<QTNode *>>(max_level+1);
    int r = num_radial_points;
    int t = num_azimuthal_points;
    for (int i = 0; i < (max_level + 1); i++) {
        levelVector[i] = std::vector<QTNode *>(r*t);
        // Divide by 2 leaving at least 1
        r = max(1, r>>1);
        t = max(1, t>>1);
    }

    // Recurseively create the node structure
    head = createNode((QTNode *) nullptr, 0, pow(2, max_level) - 1, 0, pow(2, max_level) - 1, max_level);

    // Iterate through all nodes and find their neighbours
    for (int level = 0; level < (max_level + 1); level++) {
        for (QTNode *node : levelVector[level]) {
            node->neighbour = vector<QTNode *>();

            // calculate 8 possible neigbours
            if (node->r_start > 0) {
                node->neighbour.push_back(levelVector[level][calcIndex(level, node->r_start - 1, node->t_start)] );
                node->neighbour.push_back(levelVector[level][calcIndex(level, node->r_start - 1, (node->t_end + 1) % num_azimuthal_points)] );
                node->neighbour.push_back(levelVector[level][calcIndex(level, node->r_start - 1, (node->t_start + num_azimuthal_points - 1) % num_azimuthal_points)] );
            }
            if (node->r_end < num_radial_points - 1) {
                node->neighbour.push_back(levelVector[level][calcIndex(level, node->r_end + 1, node->t_start)] );
                node->neighbour.push_back(levelVector[level][calcIndex(level, node->r_end + 1, (node->t_end + 1) % num_azimuthal_points)] );
                node->neighbour.push_back(levelVector[level][calcIndex(level, node->r_end + 1, (node->t_start + num_azimuthal_points - 1) % num_azimuthal_points)] );
            }
            node->neighbour.push_back(levelVector[level][calcIndex(level, node->r_start, (node->t_end + 1) % num_azimuthal_points)] );
            node->neighbour.push_back(levelVector[level][calcIndex(level, node->r_start, (node->t_start + num_azimuthal_points - 1) % num_azimuthal_points)] );
        }
    }
}



QTNode *QuadTree::createNode(QTNode *parent, int r_start, int r_end, int t_start, int t_end, unsigned tree_level) {

    QTNode *node = nullptr;

    // If this node is outside the range don't create it and return null
    if ((r_start < num_radial_points) &&
        (t_start < num_azimuthal_points) ) {

        unsigned step = (unsigned) pow(2, max((unsigned)0, tree_level-1));

        node = new QTNode(parent, r_start, r_end, t_start, t_end, tree_level);
        levelVector[tree_level][calcIndex(tree_level, r_start, t_start)] = node;

        // If this node is a single cell then leave the leaves null
        if (tree_level > 0) {

            node->leaf[0] = createNode(node,
                                       r_start, r_start + step - 1,
                                       t_start, t_start + step - 1,
                                       tree_level - 1);

            node->leaf[1] = createNode(node,
                                       r_start + step, r_end,
                                       t_start, t_start + step - 1,
                                       tree_level - 1);

            node->leaf[2] = createNode(node,
                                       r_start, r_start + step - 1,
                                       t_start + step, t_end,
                                       tree_level - 1);

            node->leaf[3] = createNode(node,
                                       r_start + step, r_end,
                                       t_start + step, t_end,
                                       tree_level - 1);


        }
    }
    return node;
}


void QuadTree::printOneNode(QTNode *node) {
    cout << node << "\n";
    if (node != nullptr) {
        for (int i = 0; i < 4; i++) {
            if (node->leaf[i] != nullptr) {
                printOneNode(node->leaf[i]);
            }
        }
    }
}

void QuadTree::printNodesPerLevel() {
    int level = 0;
    for (std::vector<QTNode *> nodes : levelVector) {
        cout << "Level:" << (max_level - level) << "\n";
        int linebreak = 0;
        for (QTNode *node : nodes) {
            cout << (void *) node << ",";
            if (linebreak++ > 14) {
                cout << "\n";
                linebreak = 0;
            }
        }
        cout << "\n";
        level++;
    }
}

void QuadTree::printNodes() {

    std::streambuf *psbuf, *backup;
    std::ofstream filestr;
    string filename("QuadTree.txt");
    filestr.open(filename);

    backup = std::cout.rdbuf();     // back up cout's streambuf

    psbuf = filestr.rdbuf();        // get file's streambuf
    std::cout.rdbuf(psbuf);         // assign streambuf to cout

    printOneNode(head);
    printNodesPerLevel();

    std::cout.rdbuf(backup);        // restore cout's original streambuf

    filestr.close();

    cout << "INFO: QuadTree written to " << filename;


}

QTNode *QuadTree::getHead() const {
    return head;
}

std::vector<QTNode *> QuadTree::getLevelVector(int level) {
    return levelVector[level];

}

int QuadTree::getMaxLevel() { return max_level; }