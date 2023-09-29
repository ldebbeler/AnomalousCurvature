//File: scaling.h
#ifndef SCALING_H
#define SCALING_H

#include <complex>

class scaling{
    double m_T;

public:
    scaling(double T);

private:
    std::complex<double> termx(double k0, double kt, double x);
    std::complex<double> termZero(double k0, double kt);
    std::complex<double> integrand(double k0, double kt, double x);
    std::complex<double> logIntegrand(double K, double Kt, double x);
public:
    // value of scaling function
    std::complex<double> value(double x);

    // first derivative of scaling function
    std::complex<double> firstDIntegrand(double k0, double kt, double x);
    std::complex<double> firstDlogIntegrand(double K0, double Kt, double x);
    std::complex<double> firstDvalue(double x);
    
    // second derivative of scaling function
    std::complex<double> secondDIntegrand(double k0, double kt, double x);
    std::complex<double> secondDlogIntegrand(double K0, double Kt, double x);
    std::complex<double> secondDvalue(double x);

    // coefficients of asymptotic behavior
    std::complex<double> coefficientZeroIntegrand(double k0, double kt, double sx);

    std::complex<double> coefficientZero(double sx);

    std::complex<double> coefficientTwoIntegrand(double k0, double kt, double sx);

    std::complex<double> coefficientTwo(double sx);
};
#endif // SCALING_H
