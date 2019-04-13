#include <fstream>
#include <armadillo>
#include "Grid.h"

Grid::~Grid() {
    deleteTwoDimensionalArray(nodes, columns);
    for (auto &element : elements) {
        delete element;
    }
}

Grid::Grid(std::string filename) {
    std::ifstream file(filename);
    std::stringstream buffer;
    buffer << file.rdbuf();
    file.close();

    double temperature, kFactor;
    buffer >> height >> width >> rows >> columns;
    buffer >> temperature >> kFactor;
    buffer >> heatCapacity >> density >> alpha >> ambientTemperature;

    for (auto &boundaryCondition : boundaryConditions) {
        buffer >> boundaryCondition;
    }

    init(temperature, kFactor);
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

double **Grid::calculateLocalHMatrix(const Element *element, const std::vector<JacobianMatrix> &jacobianMatrices) {
    double dndx[Constants::ELEMENT_NODE_COUNT][Constants::ELEMENT_NODE_COUNT];
    double dndy[Constants::ELEMENT_NODE_COUNT][Constants::ELEMENT_NODE_COUNT];
    //shape functions derivatives
    for (int i = 0; i < Constants::ELEMENT_NODE_COUNT; ++i) {
        double **const inversedMatrix = jacobianMatrices[i].getInversedMatrix();
        for (int j = 0; j < Constants::ELEMENT_NODE_COUNT; ++j) {
            dndx[i][j] = inversedMatrix[0][0] * Constants::dNdKsi[i][j] + inversedMatrix[0][1] * Constants::dNdEta[i][j];
            dndy[i][j] = inversedMatrix[1][0] * Constants::dNdKsi[i][j] + inversedMatrix[1][1] * Constants::dNdEta[i][j];
        }
    }

    auto H = initializeTwoDimensionalArray(Constants::ELEMENT_NODE_COUNT, Constants::ELEMENT_NODE_COUNT);
    //H matrix: integral(k * (d{N}dx * d{N}^Tdx + d{N}dy * d{N}^Tdy))dV + boundary conditions
    for (int i = 0; i < Constants::ELEMENT_NODE_COUNT; ++i) {
        double determinant = jacobianMatrices[i].getDeterminant();
        for (int j = 0; j < Constants::ELEMENT_NODE_COUNT; ++j) {
            for (int k = 0; k < Constants::ELEMENT_NODE_COUNT; ++k) {
                //values of function under integral for integral points
                double f = element->getK() * (dndx[i][j] * dndx[i][k] + dndy[i][j] * dndy[i][k]); //multiplied times wage which is 1
                //to avoid integral error, multiply function times jacobian determinant
                H[j][k] += f * determinant;
            }
        }
    }

    //boundary conditions: integral(alpha * {N} * {N}^T)dS
    auto segmentLengths = calculateSegmentLengths(element);
    auto elementEdges = element->getEdges(rows, columns);
    for (int segmentIndex = 0; segmentIndex < Constants::BOUNDARY_COUNT; ++segmentIndex) {
        if (!boundaryConditions[segmentIndex] || ((elementEdges & (1 << segmentIndex)) == 0)) {
            continue;
        }
        double determinant = segmentLengths[segmentIndex] / 2; //determinant of 1D jacobian
        for (auto &integrationPoint : Constants::SEGMENT_INTEGRATION_POINTS[segmentIndex]) {
            auto shapeFunctions = pointShapeFunctions(integrationPoint);
            for (int i = 0; i < Constants::ELEMENT_NODE_COUNT; ++i) {
                for (int j = 0; j < Constants::ELEMENT_NODE_COUNT; ++j) {
                    H[i][j] += alpha * shapeFunctions[i] * shapeFunctions[j] * determinant;
                }
            }
        }
    }
    delete[] segmentLengths;
    return H;
}

double **Grid::calculateLocalCMatrix(const std::vector<JacobianMatrix> &jacobianMatrices) {
    auto C = initializeTwoDimensionalArray(Constants::ELEMENT_NODE_COUNT, Constants::ELEMENT_NODE_COUNT);
    //C matrix: integral(c * ro * {N} * {N}^T)dV
    for (int i = 0; i < Constants::ELEMENT_NODE_COUNT; ++i) { //i - integration point index
        auto determinant = jacobianMatrices[i].getDeterminant();
        for (int j = 0; j < Constants::ELEMENT_NODE_COUNT; ++j) {
            for (int k = 0; k < Constants::ELEMENT_NODE_COUNT; ++k) {
                C[j][k] += Constants::SHAPE_FUNCTIONS[i][j] * Constants::SHAPE_FUNCTIONS[i][k] * determinant * heatCapacity * density;
            }
        }
    }
    return C;
}


double *Grid::calculateLocalPVector(const Element *element) {
    auto p = new double[Constants::ELEMENT_NODE_COUNT];
    for (int i = 0; i < Constants::ELEMENT_NODE_COUNT; ++i) {
        p[i] = 0;
    }
    //boundary conditions: integral(alpha * {N} * t)dS
    auto segmentLengths = calculateSegmentLengths(element);
    auto elementEdges = element->getEdges(rows, columns);
    for (int segmentIndex = 0; segmentIndex < Constants::BOUNDARY_COUNT; ++segmentIndex) {
        if (!boundaryConditions[segmentIndex] || ((elementEdges & (1 << segmentIndex)) == 0)) {
            continue;
        }
        double determinant = segmentLengths[segmentIndex] / 2; //determinant of 1D jacobian
        for (auto &integrationPoint : Constants::SEGMENT_INTEGRATION_POINTS[segmentIndex]) {
            auto shapeFunctions = pointShapeFunctions(integrationPoint);
            for (int i = 0; i < Constants::ELEMENT_NODE_COUNT; ++i) {
                p[i] += alpha * shapeFunctions[i] * ambientTemperature * determinant;
            }
            delete[] shapeFunctions;
        }
    }

    delete[] segmentLengths;
    return p;
}

double *Grid::calculateSegmentLengths(const Element *const element) {
    auto lengths = new double[Constants::BOUNDARY_COUNT];
    for (int i = 0; i < Constants::BOUNDARY_COUNT; ++i) {
        int j = (i + 1) % Constants::BOUNDARY_COUNT; //next segment index
        auto iNode = element->getNode(i);
        auto jNode = element->getNode(j);
        lengths[i] = sqrt(pow(jNode->getX() - iNode->getX(), 2) + pow(jNode->getY() - iNode->getY(), 2));
    }
    return lengths;
}

void Grid::evaluate(const double &time, const double &tau) {
    auto nodeCount = columns * rows;

    auto H = initializeTwoDimensionalArray(nodeCount, nodeCount);
    auto C = initializeTwoDimensionalArray(nodeCount, nodeCount);
    auto P = new double[nodeCount];
    for (int i = 0; i < nodeCount; ++i) {
        P[i] = 0;
    }

    /* Calculate parameters */
    for (auto &element : elements) {
        std::vector<JacobianMatrix> jacobianMatrices;
        jacobianMatrices.reserve(Constants::ELEMENT_NODE_COUNT);
        for (int i = 0; i < Constants::ELEMENT_NODE_COUNT; ++i) {
            jacobianMatrices.emplace_back(element, i);
        }

        auto h = calculateLocalHMatrix(element, jacobianMatrices);
        auto c = calculateLocalCMatrix(jacobianMatrices);
        auto p = calculateLocalPVector(element);

        /* Put local values to global matrices/vector */
        for (int i = 0; i < Constants::ELEMENT_NODE_COUNT; ++i) {
            auto iNodeId = element->getNode(i)->getId();
            for (int j = 0; j < Constants::ELEMENT_NODE_COUNT; ++j) {
                auto jNodeId = element->getNode(j)->getId();
                H[iNodeId][jNodeId] += h[i][j];
                C[iNodeId][jNodeId] += c[i][j];
            }
            P[iNodeId] += p[i];
        }
        deleteTwoDimensionalArray(h, Constants::ELEMENT_NODE_COUNT);
        deleteTwoDimensionalArray(c, Constants::ELEMENT_NODE_COUNT);
        delete[] p;
    }
    /* Prepare equation coefficient */
    for (int i = 0; i < nodeCount; ++i) {
        for (int j = 0; j < nodeCount; ++j) {
            C[i][j] /= tau;
        }
    }
    for (int i = 0; i < nodeCount; ++i) {
        for (int j = 0; j < nodeCount; ++j) {
            H[i][j] += C[i][j];
        }
    }

    /* Get initial temperatures */
    auto temperatures = new double[nodeCount];
    for (int i = 0; i < nodeCount; ++i) {
        temperatures[i] = nodes[i / rows][i % rows].getTemperature();
    }

    /* Simulation */
    for (int t = 0; t < ceil(time / tau); ++t) {
        std::cout << t << " iteration\n";
        auto p = new double[nodeCount];
        for (int i = 0; i < nodeCount; ++i) {
            double c = 0;
            for (int j = 0; j < nodeCount; ++j) {
                c += C[i][j] * temperatures[j];
            }
            p[i] = P[i] + c;
        }

        delete[] temperatures;
        temperatures = solveEquation(H, p);
        delete[] p;
    }

    /* Cleaning */
    deleteTwoDimensionalArray(H, nodeCount);
    deleteTwoDimensionalArray(C, nodeCount);
    delete[] P;
    delete[] temperatures;
}

template<typename T>
void Grid::deleteTwoDimensionalArray(T **array, const int &rows) {
    for (int i = 0; i < rows; ++i) {
        delete[] array[i];
    }
    delete[] array;
}

double **Grid::initializeTwoDimensionalArray(const int &n, const int &m) {
    auto array = new double *[n];
    for (int i = 0; i < n; ++i) {
        array[i] = new double[m];
        for (int j = 0; j < m; ++j) {
            array[i][j] = 0;
        }
    }
    return array;
}

double *Grid::pointShapeFunctions(const Point &point) {
    auto shapeFunctions = new double[Constants::ELEMENT_NODE_COUNT];
    shapeFunctions[0] = 0.25 * (1 - point.x) * (1 - point.y);
    shapeFunctions[1] = 0.25 * (1 + point.x) * (1 - point.y);
    shapeFunctions[2] = 0.25 * (1 + point.x) * (1 + point.y);
    shapeFunctions[3] = 0.25 * (1 - point.x) * (1 + point.y);
    return shapeFunctions;
}

double *Grid::solveEquation(double **H, double *p) {
    auto nodeCount = static_cast<unsigned long long>(columns * rows);
    arma::mat A(nodeCount, nodeCount);
    arma::vec b(nodeCount);
    for (arma::uword i = 0; i < nodeCount; ++i) {
        for (arma::uword j = 0; j < nodeCount; ++j) {
            A(i, j) = H[i][j];
        }
        b(i) = p[i];
    }
    auto resultVector = arma::solve(A, b).eval();
    std::cout << "MIN: " << resultVector.min() << "; MAX: " << resultVector.max() << std::endl;
    auto result = new double[nodeCount];
    for (unsigned long long i = 0; i < nodeCount; ++i) {
        result[i] = resultVector(i);
    }
    return result;
}