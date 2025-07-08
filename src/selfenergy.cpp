#include "selfenergy.h"
#include "constants.h"
#include "d2.hpp"
#include "integrate2d.hpp"
#include "integrategsl.hpp"
#include "linearFit.h"
#include<iostream>

selfenergy::selfenergy(double T): m_T(T), m_scaleReal(bubbleScaleReal(T)), m_scaleImag(bubbleScaleImag(T)) {}

std::complex<double> selfenergy::bubbleScale(double x){
    return m_scaleReal.evaluate(x) + I*m_scaleImag.evaluate(x);
}

std::complex<double> selfenergy::exactBubbleQuartic(double x){
    if(x<0.125){
        return 1.0/M_PI*std::sqrt(3.0/2.0+std::sqrt(std::abs(0.25-2.0*x)));
    }
    else{
        return 1.0/(2.0*M_PI)*(std::sqrt(2*std::sqrt(2.0*(1+x))+3.0) - I*std::sqrt(2.0*std::sqrt(2.0*(1+x))-3.0));
    }
}

double selfenergy::imagBubbleDeriv(double x){
    if(x<std::pow(2,1.0-m_T)){
        return 0.0;
    }
    else{
        double kzero{ -M_PI*m_scaleImag.evaluate(x) };
        auto asymm = [&] (double x) -> double{ return std::pow(std::abs(x),m_T-1.0)*sgn(x); };
        return -1.0/(M_PI*m_T*(asymm(kzero+0.5) + asymm(kzero-0.5)));
    }
}

double selfenergy::exactDerivQuartic(double x){
    if(x<1/8){
        return 0.0;
    }
    else{
        return -1.0/(2.0*M_PI*std::sqrt(2*(1+x))*std::sqrt(2*std::sqrt(2*(1+x))-3));
    }
}

std::complex<double> selfenergy::freqIntegrandPos(double qr, double qt){
    double x1{ (std::pow(std::abs(qt),m_T)-1.0)/(std::pow(std::abs(qt),m_T)) };
    double x2{ (-2*qr+std::pow(std::abs(qt),m_T)+1.0)/(std::pow(std::abs(qt),m_T)) };
    return -8.0/(std::abs(qt)*(bubbleScale(x1)+std::conj(bubbleScale(x2))));
}

std::complex<double> selfenergy::freqIntegrandNeg(double qr, double qt){
    double x1{ (std::pow(std::abs(qt),m_T)+1.0)/(std::pow(std::abs(qt),m_T)) };
    double x2{ (2*qr+std::pow(std::abs(qt),m_T)-1.0)/(std::pow(std::abs(qt),m_T)) };
    return 8.0/(std::abs(qt)*(bubbleScale(x1)+std::conj(bubbleScale(x2))));
}

double selfenergy::Apos(){
    auto integ = [&] (double2 Q) -> std::complex<double>{
        double Qr = Q.x;
        double Qt = Q.y;
        return std::exp(Qt)*std::exp(Qr)*(freqIntegrandPos(std::exp(Qr),std::exp(Qt))+freqIntegrandPos(1.0-std::exp(Qr),std::exp(Qt))).imag();
    };
    return integrate2Dmeasure2pi(std::log(IR),std::log(0.5),std::log(IR),std::log(UV)).integrate(integ).real();
}

double selfenergy::Aneg(){
    auto integ = [&] (double2 Q) -> std::complex<double>{
        double Qr = Q.x;
        double Qt = Q.y;
        return std::exp(Qt)*std::exp(Qr)*(freqIntegrandNeg(std::exp(Qr),std::exp(Qt))+freqIntegrandNeg(1.0-std::exp(Qr),std::exp(Qt))).imag();
    };
    return integrate2Dmeasure2pi(std::log(IR),std::log(0.5),std::log(IR),std::log(UV)).integrate(integ).real();
}

std::complex<double> selfenergy::radIntegrandPos(double qr, double qt, double kr){
    // factor 2 for symmetric qt integral
    double x1{ (std::pow(std::abs(qt),m_T)-kr-1.0)/(std::pow(std::abs(qt),m_T)) };
    double x2{ (-2*qr-kr+std::pow(std::abs(qt),m_T)+1.0)/(std::pow(std::abs(qt),m_T)) };
    return -8.0/(std::abs(qt)*(bubbleScale(x1)+std::conj(bubbleScale(x2))));
}

std::complex<double> selfenergy::radIntegrandNeg(double qr, double qt, double kr){
    // factor 2 for symmetric qt integral
    double x1{ (std::pow(std::abs(qt),m_T)-kr+1.0)/(std::pow(std::abs(qt),m_T)) };
    double x2{ (2*qr-kr+std::pow(std::abs(qt),m_T)-1.0)/(std::pow(std::abs(qt),m_T)) };
    return 8.0/(std::abs(qt)*(bubbleScale(x1)+std::conj(bubbleScale(x2))));
}

double selfenergy::arPos(double kr){
    auto integ = [&] (double2 Q) -> std::complex<double>{
        double Qr = Q.x;
        double Qt = Q.y;
        //return std::exp(Qt)*std::exp(Qr)*(radIntegrandPos(std::exp(Qr),std::exp(Qt),kr)+radIntegrandPos(1.0-std::exp(Qr),std::exp(Qt),kr)).imag();
        return std::exp(Qt)*(radIntegrandPos(Qr,std::exp(Qt),kr).imag());
    };
    //return integrate2Dmeasure2pi(std::log(IR),std::log(0.5),std::log(IR),std::log(UV)).integrate(integ).real();
    //return integrate2Dmeasure2pi(0.0,1.0,std::log(IR2d),std::log(UV2d)).integrate(integ).real();
    return integrate2Dmeasure2pi(0.0,1.0,std::log(IR2d),std::log(UV2d)).integrate(integ, prec2d, steps2d).real();
}

double selfenergy::arNeg(double kr){
    auto integ = [&] (double2 Q) -> std::complex<double>{
        double Qr = Q.x;
        double Qt = Q.y;
        //return std::exp(Qt)*std::exp(Qr)*(radIntegrandNeg(std::exp(Qr),std::exp(Qt),kr)+radIntegrandNeg(1.0-std::exp(Qr),std::exp(Qt),kr)).imag();
        return std::exp(Qt)*(radIntegrandNeg(Qr,std::exp(Qt),kr)).imag();
    };
    //return integrate2Dmeasure2pi(std::log(IR),std::log(0.5),std::log(IR),std::log(UV)).integrate(integ).real();
    //return integrate2Dmeasure2pi(0.0,1.0,std::log(IR2d),std::log(UV2d)).integrate(integ).real();
    return integrate2Dmeasure2pi(0.0,1.0,std::log(IR2d),std::log(UV2d)).integrate(integ, prec2d, steps2d).real();
}

interpolater1d selfenergy::radPosSpline(){
    std::vector<double> x(nodesR);
    std::vector<double> f(nodesR);
    double k{ 2*LambdaR/(nodesR-1) };
    for(int i=0; i<nodesR; i++){
        double val{ -LambdaR + k*i };
        x[i] = val;
        f[i] = arPos(val);
    }
    interpolater1d foo(x,f);
    return foo;
}

interpolater1d selfenergy::radNegSpline(){
    std::vector<double> x(nodesR);
    std::vector<double> f(nodesR);
    double k{ 2*LambdaR/(nodesR-1) };
    for(int i=0; i<nodesR; i++){
        double val{ -LambdaR + k*i };
        x[i] = val;
        f[i] = arNeg(val);
    }
    interpolater1d foo(x,f);
    return foo;
}

double selfenergy::asymptoticIntegrand(double qt, double sk){
    double x{ (std::pow(std::abs(qt),m_T)-sk)/(std::pow(std::abs(qt),m_T)) };
    return -1.0/(4*M_PI*M_PI)*imagBubbleDeriv(x)/(std::pow(std::abs(qt),m_T+1.0)*std::pow(std::abs(m_scaleReal.evaluate(x)),2.0));
}

double selfenergy::radAsymptoticPos(){
    auto integ = [&] (double Q) -> double{ return 2.0*std::exp(Q)*asymptoticIntegrand(std::exp(Q),1.0); };
    return integrater().integrate(integ, std::log(IR), std::log(UV));
}

double selfenergy::radAsymptoticNeg(){
    auto integ = [&] (double Q) -> double{ return 2.0*std::exp(Q)*asymptoticIntegrand(std::exp(Q),-1.0); };
    return integrater().integrate(integ, std::log(IR), std::log(UV));
}

std::complex<double> selfenergy::tangIntegrandPos(double qr, double qt, double kt){
    // qt Integral here not symmetric
    double x1{ (std::pow(std::abs(kt-qt),m_T)-1.0)/(std::pow(std::abs(qt),m_T)) };
    double x2{ (-2*qr+std::pow(std::abs(kt-qt),m_T)+1.0)/(std::pow(std::abs(qt),m_T)) };
    if(exactQuartic){
        return -4.0/(std::abs(qt)*(exactBubbleQuartic(x1)+std::conj(exactBubbleQuartic(x2))));
    }
    else{
        return -4.0/(std::abs(qt)*(bubbleScale(x1)+std::conj(bubbleScale(x2))));
    }
}

std::complex<double> selfenergy::tangIntegrandNeg(double qr, double qt, double kt){
    // qt Integral here not symmetric
    double x1{ (std::pow(std::abs(kt-qt),m_T)+1.0)/(std::pow(std::abs(qt),m_T)) };
    double x2{ (2*qr+std::pow(std::abs(kt-qt),m_T)-1.0)/(std::pow(std::abs(qt),m_T)) };
     if(exactQuartic){
        return 4.0/(std::abs(qt)*(exactBubbleQuartic(x1)+std::conj(exactBubbleQuartic(x2))));
    }
    else{
        return 4.0/(std::abs(qt)*(bubbleScale(x1)+std::conj(bubbleScale(x2))));
    }
}

double selfenergy::atPos(double kt){
    auto integ = [&] (double2 Q) -> std::complex<double>{
        double Qr = Q.x;
        double Qt = Q.y;
        return std::exp(Qt)*(tangIntegrandPos(Qr,std::exp(Qt),kt).imag() + tangIntegrandPos(Qr,-std::exp(Qt),kt).imag());
    };
    //return integrate2Dmeasure2pi(0.0,1.0,std::log(IR2d),std::log(UV2d)).integrate(integ).real();
    return integrate2Dmeasure2pi(0.0,1.0,std::log(IR2d),std::log(UV2d)).integrate(integ,prec2d, steps2d).real();
}

double selfenergy::atNeg(double kt){
    auto integ = [&] (double2 Q) -> std::complex<double>{
        double Qr = Q.x;
        double Qt = Q.y;
        return std::exp(Qt)*(tangIntegrandNeg(Qr,std::exp(Qt),kt).imag() + tangIntegrandNeg(Qr,-std::exp(Qt),kt).imag());
    };
    //return integrate2Dmeasure2pi(0.0,1.0,std::log(IR2d),std::log(UV2d)).integrate(integ).real();
    return integrate2Dmeasure2pi(0.0,1.0,std::log(IR2d),std::log(UV2d)).integrate(integ,prec2d, steps2d).real();
}

double selfenergy::atPosSubtract(double kt){
    auto integ = [&] (double2 Q) -> std::complex<double>{
        double Qr = Q.x;
        double Qt = Q.y;
        return std::exp(Qt)*(tangIntegrandPos(Qr, std::exp(Qt),kt).imag() - tangIntegrandPos(Qr, std::exp(Qt),0.0).imag()
                           + tangIntegrandPos(Qr,-std::exp(Qt),kt).imag() - tangIntegrandPos(Qr,-std::exp(Qt),0.0).imag());
    };
    return integrate2Dmeasure2pi(0.0,1.0,std::log(IR2d),std::log(UV2d)).integrate(integ).real();
}

double selfenergy::atNegSubtract(double kt){
    auto integ = [&] (double2 Q) -> std::complex<double>{
        double Qr = Q.x;
        double Qt = Q.y;
        return std::exp(Qt)*(tangIntegrandNeg(Qr, std::exp(Qt),kt).imag() - tangIntegrandNeg(Qr, std::exp(Qt),0.0).imag() 
                           + tangIntegrandNeg(Qr,-std::exp(Qt),kt).imag() - tangIntegrandNeg(Qr,-std::exp(Qt),0.0).imag());
    };
    return integrate2Dmeasure2pi(0.0,1.0,std::log(IR2d),std::log(UV2d)).integrate(integ).real();
}

interpolater1d selfenergy::tangPosSpline(){
    std::vector<double> x(nodesT);
    std::vector<double> f(nodesT);
    double k{ LambdaT/(nodesT-1) };
    for(int i=0; i<nodesT; i++){
        double val{ k*i };
        x[i] = val;
        f[i] = atPosSubtract(val);
    }
    interpolater1d foo(x,f);
    return foo;
}

interpolater1d selfenergy::tangNegSpline(){
    std::vector<double> x(nodesT);
    std::vector<double> f(nodesT);
    double k{ LambdaT/(nodesT-1) };
    for(int i=0; i<nodesT; i++){
        double val{ k*i };
        x[i] = val;
        f[i] = atNegSubtract(val);
    }
    interpolater1d foo(x,f);
    return foo;
}

double selfenergy::tangAsymptoticIntegrand(double qt){
    double x{ (std::pow(std::abs(1.0-qt),m_T))/(std::pow(std::abs(qt),m_T)) };
    return -1.0/(4*M_PI*M_PI)*imagBubbleDeriv(x)/(std::pow(std::abs(qt),m_T+1.0)*std::pow(std::abs(m_scaleReal.evaluate(x)),2.0));
}

// this is not subtracted yet!
double selfenergy::tangAsymptotic(){
    auto integ = [&] (double Q) -> double{ return std::exp(Q)*(tangAsymptoticIntegrand(std::exp(Q))+tangAsymptoticIntegrand(-std::exp(Q))); };
    return integrater().integrate(integ, std::log(IR), std::log(UV));
}

double selfenergy::posQuadratic(){
    // do quadratic regression for small kt
    // return coefficient. no constant term: use subtracted function
    std::vector<double> kt2(100);
    std::vector<double> atPos(100);
    for(std::size_t i=0; i<kt2.size(); i++){
        double kt{ 0.005*(i+1) };
        kt2[i] = kt*kt;
        atPos[i] = atPosSubtract(kt);
    }
    linearFit regress(kt2, atPos);
    return regress.slope_noConstant();
}

double selfenergy::negQuadratic(){
    std::vector<double> kt2(100);
    std::vector<double> atNeg(100);
    for(std::size_t i=0; i<kt2.size(); i++){
        double kt{ 0.005*(i+1) };
        kt2[i] = kt*kt;
        atNeg[i] = atNegSubtract(kt);
    }
    linearFit regress(kt2, atNeg);
    return regress.slope_noConstant();
}

double selfenergy::posQuartic(){
    // do quadratic regression for deltaA/kt^2
    // return coefficient.
    std::vector<double> kt2(quarticNodes);
    std::vector<double> atPosDivide(quarticNodes);
    for(std::size_t i=0; i<kt2.size(); i++){
        double ktsquare{ maxValue/quarticNodes*(i+1) };
        kt2[i] = ktsquare;
        atPosDivide[i] = atPosSubtract(std::sqrt(ktsquare))/(ktsquare);
    }
    linearFit regress(kt2, atPosDivide);
    return regress.slope_regular();
}

double selfenergy::negQuartic(){
    std::vector<double> kt2(quarticNodes);
    std::vector<double> atNegDivide(quarticNodes);
    for(std::size_t i=0; i<kt2.size(); i++){
        double ktsquare{ maxValue/quarticNodes*(i+1) };
        kt2[i] = ktsquare;
        atNegDivide[i] = atNegSubtract(std::sqrt(ktsquare))/(ktsquare);
    }
    linearFit regress(kt2, atNegDivide);
    return regress.slope_regular();
}

