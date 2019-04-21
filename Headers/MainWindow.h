#ifndef MES_MAINWINDOW_H
#define MES_MAINWINDOW_H

#include <QMainWindow>
#include "Solver.h"
#include "ResultsModel.h"

namespace Ui {
    class MainWindow;
}

class MainWindow : public QMainWindow {
Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = nullptr);

    ~MainWindow() override;

private slots:
    void on_simulateButton_clicked();

private:
    Ui::MainWindow *ui;
    Solver solver;
    ResultsModel model;
};

#endif //MES_MAINWINDOW_H
