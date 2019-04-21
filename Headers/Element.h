#ifndef MES_ELEMENT_H
#define MES_ELEMENT_H

#include "Node.h"
#include "Point.h"
#include "Constants.h"
#include <cmath>
#include <iostream>

class Element {
public:
    Element(Node *elementNodes[Constants::ELEMENT_NODE_COUNT], double kk);

    /**
     * Is element placed on edge of grid?
     *
     * @param gridRows
     * @param gridColumns
     * @return 4-bit mask; from right: first bit is bottom wall, second is right, third is top and fourth is left.
     */
    int getEdges(int gridRows, int gridColumns) const;


    void setId(int newValue) {
        id = newValue;
    }

    int getId() const {
        return id;
    }

    double getK() const {
        return k;
    }

    Node *getNode(int i) const {
        return nodes[i];
    }

protected:
    Node *nodes[Constants::ELEMENT_NODE_COUNT]{};
    double k{};
    int id{};
};


#endif //MES_ELEMENT_H
