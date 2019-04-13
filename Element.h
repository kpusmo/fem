#ifndef MES_ELEMENT_H
#define MES_ELEMENT_H

#include "Node.h"
#include "Point.h"
#include "Constants.h"
#include <cmath>
#include <iostream>

class Element {
public:
    Element(Node *elementNodes[Constants::ELEMENT_NODE_COUNT], double kk) : k(kk) {
        for (int i = 0; i < Constants::ELEMENT_NODE_COUNT; ++i) {
            nodes[i] = elementNodes[i];
        }
    }

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

    /**
     * Is element placed on edge of grid?
     *
     * @param gridRows
     * @param gridColumns
     * @return int 4-bit mask; from right: first bit is bottom wall, second is right, third is top and fourth is left.
     */
    int getEdges(int gridRows, int gridColumns) const {
        int result = 0;
        for (int i = 0; i < Constants::BOUNDARY_COUNT; ++i) {
            int nodeRow = nodes[i]->getRow();
            int nodeColumn = nodes[i]->getColumn();
            if (i & 1) {
                if (nodeColumn == 0 || nodeColumn == gridColumns - 1) {
                    result += (1 << i);
                }
            } else {
                if (nodeRow == 0 || nodeRow == gridRows - 1) {
                    result += (1 << i);
                }
            }
        }
        return result;
    }

protected:
    Node *nodes[Constants::ELEMENT_NODE_COUNT]{};
    double k{};
    int id{};
};


#endif //MES_ELEMENT_H
