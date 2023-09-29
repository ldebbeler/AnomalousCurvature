//File: selfenergy.h
#ifndef SELFENERGY_H
#define SELFENERGY_H

#include "bubble.h"

class selfenergy{
    double m_T;
    bubble m_bubble;
public:
    /*
    double dispersion(double qr, double qt);
    double dispersion(double kr, double qr, double kt, double qt);
    */

    selfenergy(double T, bubble foo);

    // explicite imaginary part of the self energy
    std::complex<double> integrandPos(double omega, double kr, double kt, double qr, double qt);
    std::complex<double> integrandPosLog(double omega, double kr, double kt, double qr, double Q);
    std::complex<double> integrandNeg(double omega, double kr, double kt, double qr, double qt);
    std::complex<double> integrandNegLog(double omega, double kr, double kt, double qr, double Q);
    std::complex<double> integralLog(double omega, double kr, double kt);

    double freqScalePos();
    double freqScaleNeg();

    // scalingfunctions for radial dependence
    std::complex<double> integrandPosx(double qr, double qt, double x);
    std::complex<double> integrandPosLogx(double qr, double Qt, double x);
    double integralPosx(double x);
    std::complex<double> integrandNegx(double qr, double qt, double x);
    std::complex<double> integrandNegLogx(double qr, double Qt, double x);
    double integralNegx(double x);

    // scalingfunctions for tangential dependence
    std::complex<double> integrandPosTang(double qr, double qt, double x);
    std::complex<double> integrandPosLogTang(double qr, double Qt, double x);
    double integralPosTang(double x);
    std::complex<double> integrandNegTang(double qr, double qt, double x);
    std::complex<double> integrandNegLogTang(double qr, double Qt, double x);
    double integralNegTang(double x);


};

#endif //SELFENERGY_H
