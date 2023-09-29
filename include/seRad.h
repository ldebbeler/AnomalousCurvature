//File: seRad.h
#ifndef SERAD_H
#define SERAD_H

#include "radialScale.h"

class seRad{
    radialScalePos m_scalePos;
    radialScaleNeg m_scaleNeg;

public:
    seRad(double T);
    seRad(radialScalePos scalePos, radialScaleNeg scaleNeg);

    double etaIntegrand(double kr, double omega);
    double radValue(double kr);

};

#endif //SERAD_H
