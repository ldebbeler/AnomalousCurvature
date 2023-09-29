#include "tangentialScale.h"

tangentialScalePos::tangentialScalePos(double T): newSe(T), m_T(T), m_spline(tangPosSpline()),
    m_asymptoticCoefficient(tangAsymptotic()), m_quadraticCoefficient(posQuadratic()), m_zeroValue(atPos(0.0)),
    m_quarticCoefficient(posQuartic()) {}

double tangentialScalePos::asymptoticFunction(double kt){
    return m_asymptoticCoefficient*std::pow(std::abs(kt),-m_T);
}

double tangentialScalePos::evaluate(double kt){
    double k{ std::abs(kt) };
    if(k>LambdaT-0.05){
        return asymptoticFunction(k) - m_zeroValue - m_quadraticCoefficient*k*k;
    }
    else if(k<quart){
        return m_quarticCoefficient*std::pow(k,4.0);
    }
    else{
        return m_spline.evaluate(k) - m_quadraticCoefficient*k*k;
    }
}

double tangentialScalePos::getQuartic(){
    return m_quarticCoefficient;
}

double tangentialScalePos::freqFunction(double omega){
    return evaluate(std::pow(b/omega,1.0/m_T));
}

tangentialScaleNeg::tangentialScaleNeg(double T): newSe(T), m_T(T), m_spline(tangNegSpline()),
    m_asymptoticCoefficient(tangAsymptotic()), m_quadraticCoefficient(negQuadratic()), m_zeroValue(atNeg(0.0)),
    m_quarticCoefficient(negQuartic()) {}

double tangentialScaleNeg::asymptoticFunction(double kt){
    return m_asymptoticCoefficient*std::pow(std::abs(kt),-m_T);
}

double tangentialScaleNeg::evaluate(double kt){
    double k{ std::abs(kt) };
    if(k>LambdaT-0.05){
        return asymptoticFunction(k) - m_zeroValue - m_quadraticCoefficient*k*k;
    }
    else if(k<quart){
        return m_quarticCoefficient*std::pow(k,4.0);
    }
    else{
        return m_spline.evaluate(k) - m_quadraticCoefficient*k*k;
    }
}

double tangentialScaleNeg::getQuartic(){
    return m_quarticCoefficient;
}

double tangentialScaleNeg::freqFunction(double omega){
    return evaluate(std::pow(b/omega,1.0/m_T));
}
