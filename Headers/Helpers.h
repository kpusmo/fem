#ifndef MES_HELPERS_H
#define MES_HELPERS_H

template<typename T>
void deleteTwoDimensionalArray(T **array, const int &rows) {
    for (int i = 0; i < rows; ++i) {
        delete[] array[i];
    }
    delete[] array;
}

#endif //MES_HELPERS_H
