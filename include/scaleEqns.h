//FILE: scaleEqns.h
#ifndef SCALEEQNS_H
#define SCALEEQNS_H

#include <complex>
#include "constants.h"
#include "interpolate1d.h"

class scaleEqns{
    double m_T;

public:
    scaleEqns(double T);

    double rootEqn(double kt, double x);
    double ktzero(double x);
    double imagValue(double x);
    double imagAsymptotic(double x);
    interpolater1d imagSpline();

    double realIntegrand(double kt, double x);
    double realValue(double x);
    double realAsymptotic(double x);
    interpolater1d realSpline();

};

#endif //SCALEEQNS_H
