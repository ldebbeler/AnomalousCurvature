#include "seRad.h"
#include "constants.h"
#include "integrategsl.hpp"

seRad::seRad(double T): m_scalePos(radialScalePos(T)), m_scaleNeg(radialScaleNeg(T)) {}

seRad::seRad(radialScalePos scalePos, radialScaleNeg scaleNeg): m_scalePos(scalePos), m_scaleNeg(scaleNeg) {}

double seRad::etaIntegrand(double kr, double omega){
    double x{ vF*kr/std::abs(omega) };
    return m_scaleNeg.evaluate(x)-m_scaleNeg.evaluate(0.0) - m_scalePos.evaluate(x) + m_scalePos.evaluate(0.0);
}

double seRad::radValue(double kr){
    auto integ = [&] (double w) -> double{ return M/(N*M_PI)*std::exp(w)*etaIntegrand(kr, std::exp(w)); };
    return integrater().integrate(integ, std::log(IR), std::log(UV));
}

