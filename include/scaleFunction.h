//File scaleFunction.h
#ifndef SCALEFUNCTION_H
#define SCALEFUNCTION_H

#include "interpolate1d.h"

//base class for scaling functions: real and imaginary bubble as well as scaling for self energy -> Kramers-kronig integrand

//inheriting from this class does not work as intended: segmentation faults when evaluating the spline
class scaleFunction{
    interpolater1d m_spline;

public:
    scaleFunction(interpolater1d spline);

    double evaluateSpline(double x);
    virtual double asymptoticLeft(double x) = 0;
    virtual double asymptoticRight(double x)= 0;

    double evaluate(double x);

};

#endif //SCALEFUNCTION_H
