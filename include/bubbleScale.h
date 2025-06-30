//File: bubbleScale.h
#ifndef BUBBLESCALE_H
#define BUBBLESCALE_H

//#include "scaleFunction.h"
#include "scaleEqns.h"

// these classes provide an interface to directly return real and imaginary values of the bubble scaling function
// it combines evaluation by spline with the asymptotic expressions
class bubbleScaleReal: public scaleEqns{
    interpolater1d m_spline;

public:
    bubbleScaleReal(double T);

private:
    double asymptoticLeft(double x);
    double asymptoticRight(double x);

public:
    double evaluate(double x);
};

class bubbleScaleImag: public scaleEqns{
    interpolater1d m_spline;

public:
    bubbleScaleImag(double T);

private:
    double asymptoticLeft(double x);
    double asymptoticRight(double x);

public:
    double evaluate(double x);
};

#endif //BUBBLESCALE_H
