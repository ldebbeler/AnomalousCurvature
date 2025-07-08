//File: seTang.h
#ifndef SETANG_H
#define SETANG_H

#include "tangentialScale.h"

// this class can be used to calculate the real part self energy by kramers kronig
// in practice, we do an expansion instead of the full integral
class seTang{
    double m_T;

public:
    tangentialScalePos m_scalePos;
    tangentialScaleNeg m_scaleNeg;

    double m_Apos;
    double m_Aneg;
    double m_Cpos;
    double m_Cneg;
    double m_Dpos;
    double m_Dneg;

    seTang(double T);

    double etaIntegrand(double kt, double omega);
    double tangValue(double kt);

    double bIntegrand(double omega);
    double deltab();

};

#endif //SETANG_H
