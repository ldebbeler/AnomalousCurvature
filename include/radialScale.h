//File: radialScale.h
#ifndef RADIALSCALE_H
#define RADIALSCALE_H

#include "newSe.h"

class radialScalePos: public newSe{
    interpolater1d m_spline;
    double m_rightCoefficient;
    double m_leftCoefficient;

public:

    radialScalePos(double T);

    double asymptoticLeft(double kr);
    double asymptoticRight(double kr);

    double evaluate(double kr);
};

class radialScaleNeg: public newSe{
    interpolater1d m_spline;
    double m_rightCoefficient;
    double m_leftCoefficient;

public:

    radialScaleNeg(double T);

    double asymptoticLeft(double kr);
    double asymptoticRight(double kr);

    double evaluate(double kr);
};


#endif //RADIALSCALE_H
