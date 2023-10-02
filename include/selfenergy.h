//File: newSe.h
#ifndef NEWSE_H
#define NEWSE_H

#include "bubbleScale.h"

template <typename T> int sgn(T val){
    return (T(0) < val) - (val < T(0));
}

class newSe{
public:
    double m_T;
    bubbleScaleReal m_scaleReal;
    bubbleScaleImag m_scaleImag;

    newSe(double T);
    
    std::complex<double> bubbleScale(double x);
    std::complex<double> exactBubbleQuartic(double x);
    double imagBubbleDeriv(double x);
    double exactDerivQuartic(double x);

    // frequency dependence
    std::complex<double> freqIntegrandPos(double qr, double qt);
    std::complex<double> freqIntegrandNeg(double qr, double qt);

    double Apos();
    double Aneg();

    // radial momentum dependence
    std::complex<double> radIntegrandPos(double qr, double qt, double kr);
    std::complex<double> radIntegrandNeg(double qr, double qt, double kr);

    double arPos(double kr);
    double arNeg(double kr);

    interpolater1d radPosSpline();
    interpolater1d radNegSpline();

    // same coefficient of asymptotic behavior for positive or negative frequencies
    double asymptoticIntegrand(double qt, double sk);
    double radAsymptoticPos();
    double radAsymptoticNeg();

    // tangential momentum dependence
    std::complex<double> tangIntegrandPos(double qr, double qt, double kt);
    std::complex<double> tangIntegrandNeg(double qr, double qt, double kt);

    double atPos(double kt);
    double atNeg(double kt);

    double atPosSubtract(double kt);
    double atNegSubtract(double kt);

    interpolater1d tangPosSpline();
    interpolater1d tangNegSpline();

    // asymptotic behavior same for positive and negative frequencies. Symmetric w.r.t. kt
    double tangAsymptoticIntegrand(double qt);
    double tangAsymptotic();

    double posQuadratic();
    double negQuadratic();

    double posQuartic();
    double negQuartic();
};


#endif //NEWSE_H
