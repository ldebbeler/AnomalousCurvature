//File: tangFlow.h
#ifndef TANGFLOW_H
#define TANGFLOW_H

#include "bubble.h"
#include "d2.hpp"

class tangFlow{
    bubble m_bubble;

public:
    tangFlow(bubble foo);

    std::complex<double> fraction(double lambda, double qr, double qt, int num, int denom);

    std::complex<double> fermionIntegrand(double lambda, double qr, double qt);

    std::complex<double> fullIntegrand(double lambda, double qr, double qt);

    std::complex<double> rightFlow(double lambda);

    void operator() (const std::vector<double> &x, std::vector<double>& dxdl, const double l);
};

#endif //TANGFLOW_H
