//File: seTang.h
#ifndef SETANG_H
#define SETANG_H

#include "tangentialScale.h"

class seTang{
    double m_T;

public:
    tangentialScalePos m_scalePos;
    tangentialScaleNeg m_scaleNeg;

private:
    double m_Apos;
    double m_Aneg;
    double m_Cpos;
    double m_Cneg;
    double m_Dpos;
    double m_Dneg;

public:
    seTang(double T);

    double etaIntegrand(double kt, double omega);
    double tangValue(double kt);

    double bIntegrand(double omega);
    double deltab();

};

#endif //SETANG_H
