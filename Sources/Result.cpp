#include <Result.h>

#include "Result.h"

Result::Result(const arma::vec &t) : temperatures(t) {}

const arma::vec Result::getTemperatures() const {
    return temperatures;
}
