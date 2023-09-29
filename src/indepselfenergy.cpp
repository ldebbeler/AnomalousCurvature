#include "indepselfenergy.h"
#include "grid.h"
#include "constants.h"
#include "d2.hpp"
#include "d3.hpp"
#include "integrate2d.hpp"
//#include "integrate2d_dcuhre.hpp"
#include "integrate3d.hpp"

indepselfenergy::indepselfenergy(double T): m_T(T), m_scale(scaling(m_T)),
                                        m_Cpos(m_scale.coefficientZero(1.0)), m_Cneg(m_scale.coefficientZero(-1.0)),
                                        m_Dpos(m_scale.coefficientTwo(1.0)),  m_Dneg(m_scale.coefficientTwo(-1.0)),
                                        m_x(gridVector()), m_scalingReal(realVector()), m_scalingImag(imagVector()),
                                        m_splineReal(interpolater1d(m_x,m_scalingReal)),
                                        m_splineImag(interpolater1d(m_x,m_scalingImag)) {}
                                        /*
                                        m_scalingFirstReal(realFirstVector()), m_scalingFirstImag(imagFirstVector()),
                                        m_splineFirstReal(interpolater1d(m_x,m_scalingFirstReal)), 
                                        m_splineFirstImag(interpolater1d(m_x,m_scalingFirstImag)),
                                        m_scalingSecondReal(realSecondVector()), m_scalingSecondImag(imagSecondVector()),
                                        m_splineSecondReal(interpolater1d(m_x,m_scalingSecondReal)),
                                        m_splineSecondImag(interpolater1d(m_x,m_scalingSecondImag))
                                        {}*/

std::vector<double> indepselfenergy::gridVector(){
    grid vec(minx, extraPolate, N);
    return vec.m_x;
}

std::vector<double> indepselfenergy::realVector(){
    std::vector<double> itx(m_x.size());
    for(std::size_t i=0; i<m_x.size(); i++){
        itx[i] = m_scale.value(m_x[i]).real();
    }
    return itx;
}

std::vector<double> indepselfenergy::imagVector(){
    std::vector<double> itx(m_x.size());
    for(std::size_t i=0; i<m_x.size(); i++){
        itx[i] = m_scale.value(m_x[i]).imag();
    }
    return itx;
}

/*
std::vector<double> indepselfenergy::realFirstVector(){
    std::vector<double> itx(m_x.size());
    for(std::size_t i=0; i<m_x.size(); i++){
        itx[i] = m_scale.firstDvalue(m_x[i]).real();
    }
    return itx;
}

std::vector<double> indepselfenergy::imagFirstVector(){
    std::vector<double> itx(m_x.size());
    for(std::size_t i=0; i<m_x.size(); i++){
        itx[i] = m_scale.firstDvalue(m_x[i]).imag();
    }
    return itx;
}

std::vector<double> indepselfenergy::realSecondVector(){
    std::vector<double> itx(m_x.size());
    for(std::size_t i=0; i<m_x.size(); i++){
        itx[i] = m_scale.secondDvalue(m_x[i]).real();
    }
    return itx;
}

std::vector<double> indepselfenergy::imagSecondVector(){
    std::vector<double> itx(m_x.size());
    for(std::size_t i=0; i<m_x.size(); i++){
        itx[i] = m_scale.secondDvalue(m_x[i]).imag();
    }
    return itx;
}
*/

std::complex<double> indepselfenergy::scaleFunction(double x){
    if(x<-extraPolate+0.1){
        return m_Cneg*std::pow(std::abs(x),1.0/m_T) + m_Dneg*std::pow(std::abs(x),-1.0/m_T);
    }
    else if(x > extraPolate-0.1){
        return m_Cpos*std::pow(std::abs(x),1.0/m_T) + m_Dpos*std::pow(std::abs(x),-1.0/m_T);
    }
    else{
        return m_splineReal.evaluate(x) + I*m_splineImag.evaluate(x);
    }
}

std::complex<double> indepselfenergy::scaleFirstD(double x){
    if(x<-extraPolate+0.1){
        return m_Cneg/m_T*std::pow(std::abs(x),1.0/m_T-1.0) - m_Dneg/m_T*std::pow(std::abs(x),-1.0/m_T-1.0);
    }
    else if(x > extraPolate-0.1){
        return m_Cpos/m_T*std::pow(std::abs(x),1.0/m_T-1.0) - m_Dpos/m_T*std::pow(std::abs(x),-1.0/m_T-1.0);
    }
    else{
        //return m_splineFirstReal.evaluate(x) + I*m_splineFirstImag.evaluate(x);
        return m_splineReal.firstDevaluate(x) + I*m_splineImag.firstDevaluate(x);
    }
}

std::complex<double> indepselfenergy::scaleSecondD(double x){
    if(x<-extraPolate+0.1){
        return m_Cneg/m_T*(1.0/m_T-1.0)*std::pow(std::abs(x),1.0/m_T-2.0) + m_Dneg/m_T*(1.0/m_T+1.0)*std::pow(std::abs(x),-1.0/m_T-2.0);
    }
    else if(x > extraPolate-0.1){
        return m_Cpos/m_T*(1.0/m_T-1.0)*std::pow(std::abs(x),1.0/m_T-2.0) + m_Dpos/m_T*(1.0/m_T+1.0)*std::pow(std::abs(x),-1.0/m_T-2.0);
    }
    else{
        //return m_splineSecondReal.evaluate(x) + I*m_splineSecondImag.evaluate(x);
        return m_splineReal.secondDevaluate(x) + I*m_splineImag.secondDevaluate(x);
    }
}

std::complex<double> indepselfenergy::bubble(double omega, double qr, double qt){
    return std::abs(qt)/(4.0*vF)*(scaleFunction((omega-vF*qr)/(b*std::pow(std::abs(qt),m_T)))
                     + std::conj(scaleFunction((-omega-vF*qr)/(b*std::pow(std::abs(qt),m_T)))));
}


std::complex<double> indepselfenergy::integrandPos(double omega, double kr, double kt, double qr, double qt){
    double om{ std::abs(omega) };
    return 1.0/bubble(vF*qr-om,qr-b/vF*std::pow(std::abs(qt-kt),m_T)+kr,qt);
}

std::complex<double> indepselfenergy::integrandPosLog(double omega, double kr, double kt, double qr, double Q){
    return std::exp(Q)*(integrandPos(omega,kr,kt,qr,std::exp(Q)) + integrandPos(omega,kr,kt,qr,-std::exp(Q)));
}

std::complex<double> indepselfenergy::integrandNeg(double omega, double kr, double kt, double qr, double qt){
    double om{ std::abs(omega) };
    return -1.0/bubble(vF*qr+om,qr-b/vF*std::pow(std::abs(qt-kt),m_T)+kr,qt);
}

std::complex<double> indepselfenergy::integrandNegLog(double omega, double kr, double kt, double qr, double Q){
    return std::exp(Q)*(integrandNeg(omega,kr,kt,qr,std::exp(Q)) + integrandNeg(omega,kr,kt,qr,-std::exp(Q)));
}

// no factor 2 here because integral w.r.t. qt not symmetric
std::complex<double> indepselfenergy::integralLog(double omega, double kr, double kt){
    if(omega>0.0){
        auto integ = [&] (double2 q) -> std::complex<double>{
            double qr = q.x;
            double Q = q.y;
            return integrandPosLog(omega, kr, kt, qr, Q);
        };
        return integrate2Dmeasure2pi(0.0, omega/vF, std::log(IRq), std::log(UVq)).integrate(integ, prec2dq, steps2dq);
    }
    else{
        auto integ = [&] (double2 q) -> std::complex<double>{
            double qr = q.x;
            double Q = q.y;
            return integrandNegLog(omega, kr, kt, qr, Q);
        };
        return integrate2Dmeasure2pi(omega/vF, 0.0, std::log(IRq), std::log(UVq)).integrate(integ, prec2dq, steps2dq);
    }
}

double indepselfenergy::freqScalePos(){
    return integralLog(1.0, 0.0, 0.0).imag();
}

double indepselfenergy::freqScaleNeg(){
    return integralLog(-1.0, 0.0, 0.0).imag();
}

// momentum dependence radial dependence
// global sign for definition scaling function?
std::complex<double> indepselfenergy::integrandPosx(double qr, double qt, double x){
    return -1.0/bubble(qr-1.0,1.0/vF*(qr+x-std::pow(std::abs(qt),m_T)),std::pow(b,-1.0/m_T)*qt);
}

//factor of 2 because qt integral symmetric
std::complex<double> indepselfenergy::integrandPosLogx(double qr, double Qt, double x){
    return 2.0*std::exp(Qt)*integrandPosx(qr, std::exp(Qt), x);
}

double indepselfenergy::integralPosx(double x){
    auto integ = [&] (double2 q) -> std::complex<double>{
        double qr = q.x;
        double Qt = q.y;
        return integrandPosLogx(qr, Qt, x);
    };
    return integrate2Dmeasure2pi(0.0, 1.0, std::log(IRq), std::log(UVq)).integrate(integ, prec2dq, steps2dq).imag();
}

std::complex<double> indepselfenergy::integrandNegx(double qr, double qt, double x){
    return 1.0/bubble(qr+1.0,1.0/vF*(qr+x-std::pow(std::abs(qt),m_T)),std::pow(b,-1.0/m_T)*qt);
}

std::complex<double> indepselfenergy::integrandNegLogx(double qr, double Qt, double x){
    return 2.0*std::exp(Qt)*integrandNegx(qr, std::exp(Qt), x);
}

double indepselfenergy::integralNegx(double x){
    auto integ = [&] (double2 q) -> std::complex<double>{
        double qr = q.x;
        double Qt = q.y;
        return integrandNegLogx(qr, Qt, x);
    };
    return integrate2Dmeasure2pi(-1.0, 0.0, std::log(IRq), std::log(UVq)).integrate(integ, prec2dq, steps2dq).imag();
    //return averageBZ(-1.0/vF, 0.0, std::log(IR), std::log(UV)).integrate_dcuhre(integ, prec2d);
}


// tangential momentum dependence
std::complex<double> indepselfenergy::integrandPosTang(double qr, double qt, double x){
    return -1.0/(vF*std::pow(b,1.0/m_T))*1.0/bubble(qr-1.0,1.0/vF*(qr-std::pow(std::abs(qt-x),m_T)),std::pow(b,-1.0/m_T)*qt);
}

std::complex<double> indepselfenergy::integrandPosLogTang(double qr, double Qt, double x){
    return std::exp(Qt)*(integrandPosTang(qr, std::exp(Qt), x)+integrandPosTang(qr,-std::exp(Qt), x));
}

double indepselfenergy::integralPosTang(double x){
    auto integ = [&] (double2 q) -> std::complex<double>{
        double qr = q.x;
        double Qt = q.y;
        return (integrandPosLogTang(qr, Qt, x) - integrandPosLogTang(qr, Qt, 0.0));
    };
    return integrate2Dmeasure2pi(0.0, 1.0, std::log(IRq), std::log(UVq)).integrate(integ, prec2dq, steps2dq).imag();
    //return averageBZ(0.0, 1.0/vF, std::log(IR), std::log(UV)).integrate_dcuhre(integ, prec2d).imag();
}

std::complex<double> indepselfenergy::integrandNegTang(double qr, double qt, double x){
    return 1.0/(vF*std::pow(b,1.0/m_T))*1.0/bubble(qr+1.0,1.0/vF*(qr-std::pow(std::abs(qt-x),m_T)),std::pow(b,-1.0/m_T)*qt);
}

std::complex<double> indepselfenergy::integrandNegLogTang(double qr, double Qt, double x){
    return std::exp(Qt)*(integrandNegTang(qr, std::exp(Qt), x) + integrandNegTang(qr,-std::exp(Qt), x));
}

double indepselfenergy::integralNegTang(double x){
    auto integ = [&] (double2 q) -> std::complex<double>{
        double qr = q.x;
        double Qt = q.y;
        return (integrandNegLogTang(qr, Qt, x) - integrandNegLogTang(qr, Qt, 0.0));
    };
    return integrate2Dmeasure2pi(-1.0, 0.0, std::log(IRq), std::log(UVq)).integrate(integ, prec2dq, steps2dq).imag();
    //return averageBZ(-1.0/vF, 0.0, std::log(IR), std::log(UV)).integrate_dcuhre(integ, prec2d).imag();
}

std::vector<double> indepselfenergy::getxVector(){
    return m_x;
}

std::vector<double> indepselfenergy::getReals(){
    return m_scalingReal;
}

std::vector<double> indepselfenergy::getImags(){
    return m_scalingImag;
}

/*
std::complex<double> indepselfenergy::derivPos(double qr, double qt, double kt){
    double v1{ (std::pow(std::abs(qt-kt),m_T)-1.0)/std::pow(std::abs(qt),m_T) };
    double v2{ (-2*qr+std::pow(std::abs(qt-kt),m_T)+1.0)/std::pow(std::abs(qt),m_T) };
    std::complex<double> num{ 4.0*m_T*std::pow(std::abs(qt-kt),m_T-1.0)/(std::pow(std::abs(qt),m_T))*(scaleFirstD(v1) + std::conj(scaleFirstD(v2))) };
    std::complex<double> den{ std::abs(qt)*std::pow(scaleFunction(v1) + std::conj(scaleFunction(v2)),2) };
    return num/den;
}

std::complex<double> indepselfenergy::derivPosKernel(double qr, double Qt, double kt){
    return std::exp(Qt)*(derivPos(qr,-std::exp(Qt)+kt, kt) - derivPos(qr, std::exp(Qt)+kt, kt));
}

double indepselfenergy::derivPosValue(double kt){
    auto integ = [&] (double2 q) -> std::complex<double>{
        double qr = q.x;
        double Qt = q.y;
        return derivPosKernel(qr, Qt, kt);
    };
    return integrate2Dmeasure2pi(0.0, 1.0, std::log(IRq), std::log(UVq)).integrate(integ, prec2dq, steps2dq).imag();
}

std::complex<double> indepselfenergy::derivNeg(double qr, double qt, double kt){
    double v1{ (std::pow(std::abs(qt-kt),m_T)+1.0)/std::pow(std::abs(qt),m_T) };
    double v2{ (-2*qr+std::pow(std::abs(qt-kt),m_T)-1.0)/std::pow(std::abs(qt),m_T) };
    std::complex<double> num{ -4.0*m_T*std::pow(std::abs(qt-kt),m_T-1.0)/(std::pow(std::abs(qt),m_T))*(scaleFirstD(v1) + std::conj(scaleFirstD(v2))) };
    std::complex<double> den{ std::abs(qt)*std::pow(scaleFunction(v1) + std::conj(scaleFunction(v2)),2) };
    return num/den;
}

std::complex<double> indepselfenergy::derivNegKernel(double qr, double Qt, double kt){
    return std::exp(Qt)*(derivNeg(qr,-std::exp(Qt)+kt, kt) - derivNeg(qr, std::exp(Qt)+kt, kt));
}

double indepselfenergy::derivNegValue(double kt){
    auto integ = [&] (double2 q) -> std::complex<double>{
        double qr = q.x;
        double Qt = q.y;
        return derivNegKernel(qr, Qt, kt);
    };
    return integrate2Dmeasure2pi(-1.0, 0.0, std::log(IRq), std::log(UVq)).integrate(integ, prec2dq, steps2dq).imag();
}


std::complex<double> indepselfenergy::firstTermPos(double qr, double qt, double kt){
    double v1{ (std::pow(std::abs(qt-kt),m_T)-1.0)/std::pow(std::abs(qt),m_T) };
    double v2{ (-2*qr+std::pow(std::abs(qt-kt),m_T)+1.0)/std::pow(std::abs(qt),m_T) };
    std::complex<double> num{ 4.0*m_T*(m_T-1.0)*std::pow(std::abs(qt-kt),m_T-2.0)/(std::pow(std::abs(qt),m_T))*(scaleFirstD(v1) + std::conj(scaleFirstD(v2))) };
    std::complex<double> den{ std::abs(qt)*std::pow(scaleFunction(v1) + std::conj(scaleFunction(v2)),2) };
    return num/den;
}

std::complex<double> indepselfenergy::secondTermPos(double qr, double qt, double kt){
    double v1{ (std::pow(std::abs(qt-kt),m_T)-1.0)/std::pow(std::abs(qt),m_T) };
    double v2{ (-2*qr+std::pow(std::abs(qt-kt),m_T)+1.0)/std::pow(std::abs(qt),m_T) };
    std::complex<double> num{ 4.0*m_T*m_T*std::pow(std::abs(qt-kt),2.0*m_T-2.0)/(std::pow(std::abs(qt),2.0*m_T))*(scaleSecondD(v1) + std::conj(scaleSecondD(v2))) };
    std::complex<double> den{ std::abs(qt)*std::pow(scaleFunction(v1) + std::conj(scaleFunction(v2)),2) };
    return num/den;
}

std::complex<double> indepselfenergy::thirdTermPos(double qr, double qt, double kt){
    double v1{ (std::pow(std::abs(qt-kt),m_T)-1.0)/std::pow(std::abs(qt),m_T) };
    double v2{ (-2*qr+std::pow(std::abs(qt-kt),m_T)+1.0)/std::pow(std::abs(qt),m_T) };
    std::complex<double> num{ -8.0*m_T*m_T*std::pow(std::abs(qt-kt),2.0*m_T-2.0)/(std::pow(std::abs(qt),2.0*m_T))*std::pow(scaleFirstD(v1)+std::conj(scaleFirstD(v2)),2) };
    std::complex<double> den{ std::abs(qt)*std::pow(scaleFunction(v1) + std::conj(scaleFunction(v2)),3) };
    return num/den;
}

std::complex<double> indepselfenergy::twoDintegrandPos(double qr, double qt, double kt){
    return firstTermPos(qr, qt, kt) + secondTermPos(qr, qt, kt) + thirdTermPos(qr, qt, kt);
}

std::complex<double> indepselfenergy::twoDintegrandLogPos(double qr, double Qt, double kt){
    return std::exp(Qt)*(twoDintegrandPos(qr, std::exp(Qt), kt) + twoDintegrandPos(qr,-std::exp(Qt), kt));
}

double indepselfenergy::twoDPos(double kt){
    auto integ = [&] (double2 q) -> std::complex<double>{
        double qr = q.x;
        double Qt = q.y;
        return twoDintegrandLogPos(qr, Qt, kt);
    };
    return integrate2Dmeasure2pi(0.0, 1.0, std::log(IRq), std::log(UVq)).integrate(integ, prec2dq, steps2dq).imag();
}

std::complex<double> indepselfenergy::firstTermNeg(double qr, double qt, double kt){
    double v1{ (std::pow(std::abs(qt-kt),m_T)+1.0)/std::pow(std::abs(qt),m_T) };
    double v2{ (-2*qr+std::pow(std::abs(qt-kt),m_T)-1.0)/std::pow(std::abs(qt),m_T) };
    std::complex<double> num{ -4.0*m_T*(m_T-1.0)*std::pow(std::abs(qt-kt),m_T-2.0)/(std::pow(std::abs(qt),m_T))*(scaleFirstD(v1) + std::conj(scaleFirstD(v2))) };
    std::complex<double> den{ std::abs(qt)*std::pow(scaleFunction(v1) + std::conj(scaleFunction(v2)),2) };
    return num/den;
}

std::complex<double> indepselfenergy::secondTermNeg(double qr, double qt, double kt){
    double v1{ (std::pow(std::abs(qt-kt),m_T)+1.0)/std::pow(std::abs(qt),m_T) };
    double v2{ (-2*qr+std::pow(std::abs(qt-kt),m_T)-1.0)/std::pow(std::abs(qt),m_T) };
    std::complex<double> num{ -4.0*m_T*m_T*std::pow(std::abs(qt-kt),2.0*m_T-2.0)/(std::pow(std::abs(qt),2.0*m_T))*(scaleSecondD(v1) + std::conj(scaleSecondD(v2))) };
    std::complex<double> den{ std::abs(qt)*std::pow(scaleFunction(v1) + std::conj(scaleFunction(v2)),2) };
    return num/den;
}

std::complex<double> indepselfenergy::thirdTermNeg(double qr, double qt, double kt){
    double v1{ (std::pow(std::abs(qt-kt),m_T)+1.0)/std::pow(std::abs(qt),m_T) };
    double v2{ (-2*qr+std::pow(std::abs(qt-kt),m_T)-1.0)/std::pow(std::abs(qt),m_T) };
    std::complex<double> num{ 8.0*m_T*m_T*std::pow(std::abs(qt-kt),2.0*m_T-2.0)/(std::pow(std::abs(qt),2.0*m_T))*std::pow(scaleFirstD(v1)+std::conj(scaleFirstD(v2)),2) };
    std::complex<double> den{ std::abs(qt)*std::pow(scaleFunction(v1) + std::conj(scaleFunction(v2)),3) };
    return num/den;
}

std::complex<double> indepselfenergy::twoDintegrandNeg(double qr, double qt, double kt){
    return firstTermNeg(qr, qt, kt) + secondTermNeg(qr, qt, kt) + thirdTermNeg(qr, qt, kt);
}

std::complex<double> indepselfenergy::twoDintegrandLogNeg(double qr, double Qt, double kt){
    return std::exp(Qt)*(twoDintegrandNeg(qr, std::exp(Qt), kt) + twoDintegrandNeg(qr,-std::exp(Qt), kt));
}

double indepselfenergy::twoDNeg(double kt){
    auto integ = [&] (double2 q) -> std::complex<double>{
        double qr = q.x;
        double Qt = q.y;
        return twoDintegrandLogNeg(qr, Qt, kt);
    };
    return integrate2Dmeasure2pi(-1.0, 0.0, std::log(IRq), std::log(UVq)).integrate(integ, prec2dq, steps2dq).imag();
}
std::complex<double> indepselfenergy::fourthIntegrand(double qr, double qt){
    double v1{ (std::pow(std::abs(qt),m_T)-1.0)/std::pow(std::abs(qt),m_T) };
    double v2{ (-2*qr+std::pow(std::abs(qt),m_T)+1.0)/std::pow(std::abs(qt),m_T) };
    std::complex<double> num{ scaleFirstD(v1) + std::conj(scaleFirstD(v2)) };
    std::complex<double> den{ std::pow(std::abs(qt),5.0)*std::pow(scaleFunction(v1) + std::conj(scaleFunction(v2)),2) };
    return num/den;
}

double indepselfenergy::fourthDerivative(double irLambda, double uvLambda){
    auto integ = [&] (double2 q) -> std::complex<double>{
        double qr = q.x;
        double Qt = q.y;
        return std::exp(Qt)*fourthIntegrand(qr,std::exp(Qt));
    };
    return integrate2Dmeasure2pi(0.0, 1.0, std::log(irLambda), std::log(uvLambda)).integrate(integ, prec2dq, steps2dq).imag();
}

std::vector<double> indepselfenergy::getFirstDReal(){
    return m_scalingFirstReal;
}

std::vector<double> indepselfenergy::getFirstDImag(){
    return m_scalingFirstImag;
}

std::vector<double> indepselfenergy::getSecondDReal(){
    return m_scalingSecondReal;
}

std::vector<double> indepselfenergy::getSecondDImag(){
    return m_scalingSecondImag;
}
*/

