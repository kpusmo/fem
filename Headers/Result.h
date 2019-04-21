#ifndef MES_RESULT_H
#define MES_RESULT_H

#include <armadillo>

class Result {
public:
    explicit Result(const arma::vec &temperatures);

    const arma::vec getTemperatures() const;

protected:
    arma::vec temperatures;
};


#endif //MES_RESULT_H
