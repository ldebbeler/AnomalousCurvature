//File: tangentialScale.h
#ifndef TANGENTIALSCALE_H
#define TANGENTIALSCALE_H

#include "selfenergy.h"

class tangentialScalePos: public selfenergy{
    double m_T;
public:
    double m_quart;
private:
    interpolater1d m_spline;        //spline through (0,0)
    double m_asymptoticCoefficient;
    double m_quadraticCoefficient;
    double m_zeroValue;
    double m_quarticCoefficient;

public:

    tangentialScalePos(double T);

    double asymptoticFunction(double kt);

    double evaluate(double kt);

    double getQuartic();
    double getQuadratic();

    double freqFunction(double omega);
};

class tangentialScaleNeg: public selfenergy{
    double m_T;
    double m_quart;
    interpolater1d m_spline;        // spline through (0,0)
    double m_asymptoticCoefficient;
    double m_quadraticCoefficient;
    double m_zeroValue;
    double m_quarticCoefficient;

public:

    tangentialScaleNeg(double T);

    double asymptoticFunction(double kt);

    double evaluate(double kt);

    double getQuartic();
    double getQuadratic();

    double freqFunction(double omega);
};

#endif //TANGENTIALSCALE_H
