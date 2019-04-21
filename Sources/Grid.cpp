#include <fstream>
#include <Grid.h>
#include <Helpers.h>

Grid::Grid(double w, double h, int c, int r) : width(w), height(h), columns(c), rows(r) {}

Grid::~Grid() {
    deleteTwoDimensionalArray(nodes, rows);
    for (auto &element : elements) {
        delete element;
    }
}

void Grid::init(double elementInitialTemperature, double elementKFactor) {
    nodes = new Node *[rows];
    double deltaH = height / (rows - 1);
    double deltaW = width / (columns - 1);
    for (int i = rows - 1; i >= 0; --i) {
        nodes[i] = new Node[columns];
        for (int j = 0; j < columns; ++j) {
            nodes[i][j].setX(j * deltaW);
            nodes[i][j].setY(i * deltaH);
            nodes[i][j].setColumn(j);
            nodes[i][j].setRow(i);
            nodes[i][j].setTemperature(elementInitialTemperature);
            if (i != rows - 1 && j != columns - 1) {
                Node *elementNodes[] = {&nodes[i][j], &nodes[i][j + 1], &nodes[i + 1][j + 1], &nodes[i + 1][j]};
                elements.insert(elements.begin() + j, new Element(elementNodes, elementKFactor));
            }
        }
    }
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < columns; ++j) {
            nodes[i][j].setId(j * rows + i);
        }
    }
    for (int i = 0; i < rows - 1; ++i) {
        for (int j = 0; j < columns - 1; ++j) {
            elements[i * (columns - 1) + j]->setId(j * (rows - 1) + i);
        }
    }
}

void Grid::drawGrid(std::ostream &out) {
    out << "Node temperatures:\n\n";
    for (int i = rows - 1; i >= 0; --i) {
        for (int j = 0; j < columns - 1; ++j) {
            out << nodes[i][j].getTemperature() << "----";
        }
        out << nodes[i][columns - 1].getTemperature() << "\n\n";
    }
    out << "Node ids:\n\n";
    for (int i = rows - 1; i >= 0; --i) {
        for (int j = 0; j < columns - 1; ++j) {
            out << nodes[i][j].getId() << "----";
        }
        out << nodes[i][columns - 1].getId() << "\n\n";
    }
    out << "Node coords:\n\n";
    for (int i = rows - 1; i >= 0; --i) {
        for (int j = 0; j < columns - 1; ++j) {
            out << nodes[i][j].getX() << " " << nodes[i][j].getY() << "----";
        }
        out << nodes[i][columns - 1].getX() << " " << nodes[i][columns - 1].getY() << "\n\n";
    }
    out << "Element ks:\n\n";
    for (int i = rows - 2; i >= 0; --i) {
        for (int j = 0; j < columns - 2; ++j) {
            out << elements[i * (columns - 1) + j]->getK() << "----";
        }
        out << elements[i * (columns - 1) + columns - 2]->getK() << "\n\n";
    }
    out << "Element node ids:\n\n";
    for (int i = rows - 2; i >= 0; --i) {
        for (int j = 0; j < columns - 2; ++j) {
            out << elements[i * (columns - 1) + j]->getNode(3)->getId() << ' ' << elements[i * (columns - 1) + j]->getNode(2)->getId() << "----";
        }
        out << elements[i * (columns - 1) + columns - 2]->getNode(3)->getId() << ' ' << elements[i * (columns - 1) + columns - 2]->getNode(2)->getId() << "\n\n";

        for (int j = 0; j < columns - 2; ++j) {
            out << elements[i * (columns - 1) + j]->getNode(0)->getId() << ' ' << elements[i * (columns - 1) + j]->getNode(1)->getId() << "----";
        }
        out << elements[i * (columns - 1) + columns - 2]->getNode(0)->getId() << ' ' << elements[i * (columns - 1) + columns - 2]->getNode(1)->getId() << "\n\n\n";
    }
    out << "Element ids:\n\n";
    for (int i = rows - 2; i >= 0; --i) {
        for (int j = 0; j < columns - 2; ++j) {
            out << elements[i * (columns - 1) + j]->getId() << "----";
        }
        out << elements[i * (columns - 1) + columns - 2]->getId() << "\n\n";
    }
    out << "Element edges:\n\n";
    for (int i = rows - 2; i >= 0; --i) {
        for (int j = 0; j < columns - 1; ++j) {
            out << elements[i * (columns - 1) + j]->getEdges(rows, columns) << ' ';
        }
        out << "\n\n";
    }
    out << "Height: " << height << "; Width: " << width << std::endl;
}

int Grid::getRows() const {
    return rows;
}

int Grid::getColumns() const {
    return columns;
}

std::vector<Element *> Grid::getElements() const {
    return elements;
}

Node **Grid::getNodes() const {
    return nodes;
}
