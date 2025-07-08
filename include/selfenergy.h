//File: selfenergy.h
#ifndef SELFENERGY_H
#define SELFENERGY_H

#include "bubbleScale.h"

template <typename T> int sgn(T val){
    return (T(0) < val) - (val < T(0));
}

class selfenergy{
    double m_T;
    bubbleScaleReal m_scaleReal;
    bubbleScaleImag m_scaleImag;

public:
    selfenergy(double T);
    
    // return complex value of bubble scaling function
    std::complex<double> bubbleScale(double x);
    // returns exact form of known function for alpha=4
    std::complex<double> exactBubbleQuartic(double x);
    // return first derivative of imaginary part, required for asymptotic formulas
private:
    double imagBubbleDeriv(double x);
    double exactDerivQuartic(double x);

    // frequency dependence
    // integrands for coefficient
    std::complex<double> freqIntegrandPos(double qr, double qt);
    std::complex<double> freqIntegrandNeg(double qr, double qt);

public:
    // frequency coefficients
    double Apos();
    double Aneg();

private:
    // radial momentum dependence
    // integrands to calculate radial scaling function
    std::complex<double> radIntegrandPos(double qr, double qt, double kr);
    std::complex<double> radIntegrandNeg(double qr, double qt, double kr);

    // radial scaling function for self energy
public:
    double arPos(double kr);
    double arNeg(double kr);

    // spline objects for positive or negative frequencies and radial momentum dependence
    interpolater1d radPosSpline();
    interpolater1d radNegSpline();

    // same coefficient of asymptotic behavior for positive or negative frequencies
    // required for full integral instead of expansion
    double asymptoticIntegrand(double qt, double sk);
    double radAsymptoticPos();
    double radAsymptoticNeg();

private:
    // tangential momentum dependence
    std::complex<double> tangIntegrandPos(double qr, double qt, double kt);
    std::complex<double> tangIntegrandNeg(double qr, double qt, double kt);

public:
    // tangential momentum scaling function for self energy
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


#endif //SELFENERGY_H
