//FILE: scaleEqns.h
#ifndef SCALEEQNS_H
#define SCALEEQNS_H

#include <complex>
#include "constants.h"
#include "interpolate1d.h"

// this class calculates the scaling function for the ph bubble
class scaleEqns{
    double m_T;

public:
    scaleEqns(double T);

private:
    double rootEqn(double kt, double x);
    double ktzero(double x);
    // direct imaginary value for argument x
    double imagValue(double x);
public:
    // asymptotic behavior for large absolute values
    double imagAsymptotic(double x);
    // spline object for intermediate values
    interpolater1d imagSpline();

private:
    double realIntegrand(double kt, double x);
    // direct real value for argument x
    double realValue(double x);
public:
    // asymptotic behavior for large absolute values
    double realAsymptotic(double x);
    // spline object for intermediate values
    interpolater1d realSpline();

};

#endif //SCALEEQNS_H
