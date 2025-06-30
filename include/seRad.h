//File: seRad.h
#ifndef SERAD_H
#define SERAD_H

#include "radialScale.h"

// this class allows for direct calculation of the real part self energy by evaluating the Kramers-Kronig integral
// in practice, expanding the radial scaling function for small arguments and integrating by hand is easier,
// but this class can be used to verify that the results agree
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
