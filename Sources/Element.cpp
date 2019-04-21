#include "Element.h"

Element::Element(Node **elementNodes, double kk) : k(kk) {
    for (int i = 0; i < Constants::ELEMENT_NODE_COUNT; ++i) {
        nodes[i] = elementNodes[i];
    }
}

int Element::getEdges(int gridRows, int gridColumns) const {
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
