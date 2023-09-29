//File: seTang.h
#ifndef SETANG_H
#define SETANG_H

#include "tangentialScale.h"

class seTang{
public:
    double m_T;
    tangentialScalePos m_scalePos;
    tangentialScaleNeg m_scaleNeg;

    seTang(double T);

    double etaIntegrand(double kt, double omega);
    double tangValue(double kt);

};

#endif //SETANG_H
