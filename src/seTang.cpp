#include "seTang.h"
#include "constants.h"
#include "integrategsl.hpp"

seTang::seTang(double T): m_T(T), m_scalePos(tangentialScalePos(T)), m_scaleNeg(tangentialScaleNeg(T)) {}

double seTang::etaIntegrand(double kt, double omega){
    double x{ std::pow(std::abs(b/omega),1.0/m_T)*kt };
    return m_scaleNeg.evaluate(x) - m_scalePos.evaluate(x);
}

double seTang::tangValue(double kt){
    auto integ = [&] (double w) -> double{ return M/(N*M_PI)*std::exp(w)*etaIntegrand(kt, std::exp(w)); };
    return integrater().integrate(integ, std::log(IRTang), std::log(UVTang));
}
