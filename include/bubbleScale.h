//File: bubbleScale.h
#ifndef BUBBLESCALE_H
#define BUBBLESCALE_H

#include "scaleFunction.h"
#include "scaleEqns.h"

class bubbleScaleReal: public scaleEqns{
    interpolater1d m_spline;

public:
    bubbleScaleReal(double T);

    double asymptoticLeft(double x);
    double asymptoticRight(double x);

    double evaluate(double x);
};

class bubbleScaleImag: public scaleEqns{
    interpolater1d m_spline;

public:
    bubbleScaleImag(double T);

    double asymptoticLeft(double x);
    double asymptoticRight(double x);

    double evaluate(double x);
};

#endif //BUBBLESCALE_H
