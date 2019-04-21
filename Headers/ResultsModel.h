#ifndef MES_RESULTSMODEL_H
#define MES_RESULTSMODEL_H

#include <QAbstractTableModel>
#include <vector>
#include "Result.h"

class ResultsModel : public QAbstractTableModel {
public:
    explicit ResultsModel(QObject *parent = nullptr);

    int rowCount(const QModelIndex &index = QModelIndex()) const override;

    int columnCount(const QModelIndex &index = QModelIndex()) const override;

    QVariant data(const QModelIndex &index, int role = Qt::DisplayRole) const override;

    void setResults(const std::vector<Result> &resultsVector, unsigned rows, unsigned columns);

    void nextResult();

    void previousResult();

protected:
    std::vector<Result> results;
    unsigned rows{};
    unsigned columns{};
    unsigned currentResult{};
};


#endif //MES_RESULTSMODEL_H
