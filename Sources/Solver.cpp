#include "Solver.h"
#include <Helpers.h>


Solver::~Solver() {
    delete grid;
}

void Solver::setParameters(double heatCapacity, double density, double alpha, double ambientTemperature) {
    this->heatCapacity = heatCapacity;
    this->density = density;
    this->alpha = alpha;
    this->ambientTemperature = ambientTemperature;
}

void Solver::setBoundaryConditions(bool bottom, bool right, bool top, bool left) {
    boundaryConditions[0] = bottom;
    boundaryConditions[1] = right;
    boundaryConditions[2] = top;
    boundaryConditions[3] = left;
}

void Solver::setGrid(Grid *g) {
    delete grid;
    grid = g;
}

double **Solver::calculateLocalHMatrix(const Element *element, const std::vector<JacobianMatrix> &jacobianMatrices) {
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
    auto elementEdges = element->getEdges(grid->getRows(), grid->getColumns());
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

double **Solver::calculateLocalCMatrix(const std::vector<JacobianMatrix> &jacobianMatrices) {
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


double *Solver::calculateLocalPVector(const Element *element) {
    auto p = new double[Constants::ELEMENT_NODE_COUNT];
    for (int i = 0; i < Constants::ELEMENT_NODE_COUNT; ++i) {
        p[i] = 0;
    }
    //boundary conditions: integral(alpha * {N} * t)dS
    auto segmentLengths = calculateSegmentLengths(element);
    auto elementEdges = element->getEdges(grid->getRows(), grid->getColumns());
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

double *Solver::calculateSegmentLengths(const Element *const element) {
    auto lengths = new double[Constants::BOUNDARY_COUNT];
    for (int i = 0; i < Constants::BOUNDARY_COUNT; ++i) {
        int j = (i + 1) % Constants::BOUNDARY_COUNT; //next segment index
        auto iNode = element->getNode(i);
        auto jNode = element->getNode(j);
        lengths[i] = sqrt(pow(jNode->getX() - iNode->getX(), 2) + pow(jNode->getY() - iNode->getY(), 2));
    }
    return lengths;
}

std::vector<Result> Solver::evaluate(const double &time, const double &tau) {
    auto nodeCount = grid->getColumns() * grid->getRows();

    auto H = initializeTwoDimensionalArray(nodeCount, nodeCount);
    auto C = initializeTwoDimensionalArray(nodeCount, nodeCount);
    auto P = new double[nodeCount];
    for (int i = 0; i < nodeCount; ++i) {
        P[i] = 0;
    }

    /* Calculate parameters */
    for (auto &element : grid->getElements()) {
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
        temperatures[i] = grid->getNodes()[i / grid->getRows()][i % grid->getRows()].getTemperature();
    }

    std::vector<Result> results;
    results.emplace_back(arma::vec(temperatures, static_cast<const arma::uword>(nodeCount)));
    delete[] temperatures;

    /* Simulation */
    for (int t = 0; t < ceil(time / tau); ++t) {
        auto p = new double[nodeCount];
        for (int i = 0; i < nodeCount; ++i) {
            double c = 0;
            for (int j = 0; j < nodeCount; ++j) {
                c += C[i][j] * results[t].getTemperatures()(j);
            }
            p[i] = P[i] + c;
        }

        results.emplace_back(solveEquation(H, p));
        delete[] p;
    }

    /* Cleaning */
    deleteTwoDimensionalArray(H, nodeCount);
    deleteTwoDimensionalArray(C, nodeCount);
    delete[] P;
    auto it = results.begin();
    return results;
}

double **Solver::initializeTwoDimensionalArray(const int &n, const int &m) {
    auto array = new double *[n];
    for (int i = 0; i < n; ++i) {
        array[i] = new double[m];
        for (int j = 0; j < m; ++j) {
            array[i][j] = 0;
        }
    }
    return array;
}

double *Solver::pointShapeFunctions(const Point &point) {
    auto shapeFunctions = new double[Constants::ELEMENT_NODE_COUNT];
    shapeFunctions[0] = 0.25 * (1 - point.x) * (1 - point.y);
    shapeFunctions[1] = 0.25 * (1 + point.x) * (1 - point.y);
    shapeFunctions[2] = 0.25 * (1 + point.x) * (1 + point.y);
    shapeFunctions[3] = 0.25 * (1 - point.x) * (1 + point.y);
    return shapeFunctions;
}

arma::vec Solver::solveEquation(double **H, const double *p) {
    auto nodeCount = static_cast<unsigned long long>(grid->getColumns() * grid->getRows());
    arma::mat A(nodeCount, nodeCount);
    arma::vec b(nodeCount);
    for (arma::uword i = 0; i < nodeCount; ++i) {
        for (arma::uword j = 0; j < nodeCount; ++j) {
            A(i, j) = H[i][j];
        }
        b(i) = p[i];
    }
    auto resultVector = arma::solve(A, b).eval();
    return resultVector;
}
