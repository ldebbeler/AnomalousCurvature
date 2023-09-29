//File: linearFit.h
#ifndef LINEARFIT_H
#define LINEARFIT_H

#include <vector>
#include <gsl/gsl_fit.h>

class linearFit{
    std::vector<double> m_x;
    std::vector<double> m_y;

public:
    linearFit(std::vector<double> x, std::vector<double> y);

    double slope_noConstant();

    double slope_regular();
    double intercept_regular();

};

#endif //LINEARFIT_H
