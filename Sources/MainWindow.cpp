#include <MainWindow.h>
#include "ui_MainWindow.h"

MainWindow::MainWindow(QWidget *parent) : QMainWindow(parent), ui(new Ui::MainWindow) {
    ui->setupUi(this);
    ui->resultsView->setModel(&model);
}

MainWindow::~MainWindow() {
    delete ui;
}

void MainWindow::on_simulateButton_clicked() {
    auto width = ui->widthInput->value();
    auto height = ui->heightInput->value();
    auto columns = static_cast<unsigned>(ui->horizontalNodeCountInput->value());
    auto rows = static_cast<unsigned>(ui->verticalNodeCountInput->value());
    Grid *grid = new Grid(width, height, columns, rows);

    auto initialTemperature = ui->initialTemperaturInput->value();
    auto kFactor = ui->kFactorInput->value();
    grid->init(initialTemperature, kFactor);
    solver.setGrid(grid);

    auto heatCapacity = ui->heapCapacityInput->value();
    auto density = ui->densityInput->value();
    auto alpha = ui->alphaInput->value();
    auto ambientTemperature = ui->ambientTemperatureInput->value();
    solver.setParameters(heatCapacity, density, alpha, ambientTemperature);

    auto bcTop = ui->topBC->isChecked();
    auto bcRight = ui->rightBC->isChecked();
    auto bcBottom = ui->bottomBC->isChecked();
    auto bcLeft = ui->leftBC->isChecked();
    solver.setBoundaryConditions(bcBottom, bcRight, bcTop, bcLeft);

    auto time = ui->timeInput->value();
    auto timeStep = ui->timeStepInput->value();

    auto results = solver.evaluate(time, timeStep);
    model.setResults(results, rows, columns);
    ui->resultsView->resizeColumnsToContents();
}
