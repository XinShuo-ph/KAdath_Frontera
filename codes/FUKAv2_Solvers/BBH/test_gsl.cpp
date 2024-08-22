#include "gsl/gsl_permutation.h"
#include <iostream>

int main() {
    gsl_permutation * p = gsl_permutation_alloc(10);
    gsl_permutation_free(p);
    std::cout << "GSL is working!" << std::endl;
    return 0;
}
