#include "MainWindow.h"
#include <QApplication>

int main(int argc, char *argv[]) {
//    Grid grid("indata/test");
//    grid.drawGrid(std::cout);
//    grid.evaluate(500, 50);
    QApplication a(argc, argv);
    MainWindow w;
    w.show();

    return a.exec();
}