#include "seTang.h"
#include "constants.h"
#include "integrategsl.hpp"

seTang::seTang(double T): m_T(T), m_scalePos(tangentialScalePos(T)), m_scaleNeg(tangentialScaleNeg(T)),
    m_Apos(m_scalePos.Apos()), m_Aneg(m_scaleNeg.Aneg()), m_Cpos(m_scalePos.getQuadratic()),
    m_Cneg(m_scaleNeg.getQuadratic()), m_Dpos(m_scalePos.getQuartic()), m_Dneg(m_scaleNeg.getQuartic()) {}

double seTang::etaIntegrand(double kt, double omega){
    double x{ std::pow(std::abs(b/omega),1.0/m_T)*kt };
    return m_scaleNeg.evaluate(x) - m_scalePos.evaluate(x);
}

double seTang::tangValue(double kt){
    auto integ = [&] (double w) -> double{ return M/(N*M_PI)*std::exp(w)*etaIntegrand(kt, std::exp(w)); };
    return integrater().integrate(integ, std::log(IRTang), std::log(UVTang));
}

double seTang::bIntegrand(double omega){
    double x{ std::pow(omega,-1.0/m_T) };
    return m_scaleNeg.evaluate(x) - m_scalePos.evaluate(x);
}

double seTang::deltab(){
    if(m_T > 3.99 || m_T < 2.01){ return 0.0; }
    double low{ std::pow(LambdaT, -m_T) };
    double high{ std::pow(m_scalePos.m_quart,-m_T) };
    double uvContribution{ M/(N*M_PI)*(m_Dneg-m_Dpos)/(4.0/m_T-1.0)*std::pow(high,1.0-4.0/m_T) };
    double irContribution{ M/(N*M_PI)*((m_Apos-m_Aneg)*low + (m_Cpos-m_Cneg)/(1.0-2.0/m_T)*std::pow(low,1.0-2.0/m_T)) };
    auto integ = [&] (double omega) -> double{ return M/(N*M_PI)*bIntegrand(omega); };
    double intermediateContribution{ integrater().integrate(integ, low, high) };
    return uvContribution + irContribution + intermediateContribution;
}
