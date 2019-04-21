
#include <ResultsModel.h>

#include "ResultsModel.h"

ResultsModel::ResultsModel(QObject *parent) : QAbstractTableModel(parent) {}

int ResultsModel::rowCount(const QModelIndex &index) const {
    return rows;
}

int ResultsModel::columnCount(const QModelIndex &index) const {
    return columns;
}

QVariant ResultsModel::data(const QModelIndex &index, int role) const {
    auto row = static_cast<unsigned int>(index.row());
    auto column = static_cast<unsigned int>(index.column());

//    if (role == Qt::BackgroundRole) {
//        return grid[row][column].getColor();
//    }

    if (role == Qt::DisplayRole) {
        return results[currentResult].getTemperatures()(row * columns + column);
    }
    return QVariant();
}

void ResultsModel::nextResult() {
    if (currentResult + 1 == results.size()) {
        return;
    }
    ++currentResult;
}

void ResultsModel::previousResult() {
    if (currentResult == 0) {
        return;
    }
    --currentResult;
}

void ResultsModel::setResults(const std::vector<Result> &resultsVector, unsigned r, unsigned c) {
    beginResetModel();
    rows = r;
    columns = c;
    results = resultsVector;
    currentResult = results.size() - 1;
    endResetModel();
}
