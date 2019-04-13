#ifndef MES_NODE_H
#define MES_NODE_H

class Node {
public:
    Node() = default;

    Node(double xx, double yy, double t) : x(xx), y(yy), temperature(t) {}

    double getX() const {
        return x;
    }

    int getId() const {
        return id;
    }

    void setId(int newValue) {
        id = newValue;
    }

    void setX(double newValue) {
        x = newValue;
    }

    double getY() const {
        return y;
    }

    void setY(double newValue) {
        y = newValue;
    }

    int getColumn() const {
        return column;
    }

    void setColumn(int newValue) {
        column = newValue;
    }

    int getRow() const {
        return row;
    }

    void setRow(int newValue) {
        row = newValue;
    }

    double getTemperature() const {
        return temperature;
    }

    void setTemperature(double newValue) {
        temperature = newValue;
    }

protected:
    double x{};
    double y{};
    double temperature{};
    int id{};
    int column{};
    int row{};
};


#endif //MES_NODE_H
