#include "radialScale.h"

radialScalePos::radialScalePos(double T): selfenergy(T), m_spline(radPosSpline()), 
    m_rightCoefficient(radAsymptoticPos()), m_leftCoefficient(radAsymptoticNeg()) {}

double radialScalePos::asymptoticLeft(double kr){
    return m_leftCoefficient*std::pow(std::abs(kr),-1.0);
}

double radialScalePos::asymptoticRight(double kr){
    return m_rightCoefficient*std::pow(std::abs(kr),-1.0);
}

double radialScalePos::evaluate(double kr){
    if(kr<-LambdaR+0.1){
        return asymptoticLeft(kr);
    }
    else if(kr>LambdaR-0.1){
        return asymptoticRight(kr);
    }
    else{
        return m_spline.evaluate(kr);
    }
}

radialScaleNeg::radialScaleNeg(double T): selfenergy(T), m_spline(radNegSpline()), 
    m_rightCoefficient(radAsymptoticPos()), m_leftCoefficient(radAsymptoticNeg()) {}

double radialScaleNeg::asymptoticLeft(double kr){
    return m_leftCoefficient*std::pow(std::abs(kr),-1.0);
}

double radialScaleNeg::asymptoticRight(double kr){
    return m_rightCoefficient*std::pow(std::abs(kr),-1.0);
}

double radialScaleNeg::evaluate(double kr){
    if(kr<-LambdaR+0.1){
        return asymptoticLeft(kr);
    }
    else if(kr>LambdaR-0.1){
        return asymptoticRight(kr);
    }
    else{
        return m_spline.evaluate(kr);
    }
}
