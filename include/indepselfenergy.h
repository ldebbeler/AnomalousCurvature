//File: indepSelfEnergy.h
#ifndef INDEPSELFENERGY_H
#define INDEPSELFENERGY_H

#include "scaling.h"
#include<vector>
#include "interpolate1d.h"

class indepselfenergy{
    double m_T;
    scaling m_scale;
    std::complex<double> m_Cpos;
    std::complex<double> m_Cneg;
    std::complex<double> m_Dpos;
    std::complex<double> m_Dneg;

    std::vector<double> m_x;
    std::vector<double> m_scalingReal;
    std::vector<double> m_scalingImag;
    interpolater1d m_splineReal;
    interpolater1d m_splineImag;

    /*
    std::vector<double> m_scalingFirstReal;
    std::vector<double> m_scalingFirstImag;
    interpolater1d m_splineFirstReal;
    interpolater1d m_splineFirstImag;

    std::vector<double> m_scalingSecondReal;
    std::vector<double> m_scalingSecondImag;
    interpolater1d m_splineSecondReal;
    interpolater1d m_splineSecondImag;
    */

    std::vector<double> gridVector();
    std::vector<double> realVector();
    std::vector<double> imagVector();

    /*
    std::vector<double> realFirstVector();
    std::vector<double> imagFirstVector();

    std::vector<double> realSecondVector();
    std::vector<double> imagSecondVector();
    */

public:
    indepselfenergy(double T);

    std::complex<double> scaleFunction(double x);
    std::complex<double> scaleFirstD(double x);
    std::complex<double> scaleSecondD(double x);

    std::complex<double> bubble(double omega, double qr, double qt);
    
    // explicite imaginary part of the self energy
    std::complex<double> integrandPos(double omega, double kr, double kt, double qr, double qt);
    std::complex<double> integrandPosLog(double omega, double kr, double kt, double qr, double Q);
    std::complex<double> integrandNeg(double omega, double kr, double kt, double qr, double qt);
    std::complex<double> integrandNegLog(double omega, double kr, double kt, double qr, double Q);
    std::complex<double> integralLog(double omega, double kr, double kt);

    // anomalous dimensions frequency
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
    double integralPosTang(double x);   // important
    std::complex<double> integrandNegTang(double qr, double qt, double x);
    std::complex<double> integrandNegLogTang(double qr, double Qt, double x);
    double integralNegTang(double x);   // important

    std::vector<double> getxVector();
    std::vector<double> getReals();
    std::vector<double> getImags();


    // first derivative of tangential scaling function
    /*
    std::complex<double> derivPos(double qr, double qt, double kt);
    std::complex<double> derivPosKernel(double qr, double Qt, double kt);
    double derivPosValue(double kt);

    std::complex<double> derivNeg(double qr, double qt, double kt);
    std::complex<double> derivNegKernel(double qr, double Qt, double kt);
    double derivNegValue(double kt);

    // second derivative of tangential scaling function
    std::complex<double> firstTermPos(double qr, double qt, double kt);
    std::complex<double> secondTermPos(double qr, double qt, double kt);
    std::complex<double> thirdTermPos(double qr, double qt, double kt);
    std::complex<double> twoDintegrandPos(double qr, double qt, double kt);
    std::complex<double> twoDintegrandLogPos(double Qr, double Qt, double kt);
    double twoDPos(double kt);

    std::complex<double> firstTermNeg(double qr, double qt, double kt);
    std::complex<double> secondTermNeg(double qr, double qt, double kt);
    std::complex<double> thirdTermNeg(double qr, double qt, double kt);
    std::complex<double> twoDintegrandNeg(double qr, double qt, double kt);
    std::complex<double> twoDintegrandLogNeg(double Qr, double Qt, double kt);
    double twoDNeg(double kt);

    std::complex<double> fourthIntegrand(double qr, double qt);
    double fourthDerivative(double irlambda, double uvlambda);
    std::vector<double> getFirstDReal();
    std::vector<double> getFirstDImag();
    std::vector<double> getSecondDReal();
    std::vector<double> getSecondDImag();
*/

};

#endif //INDEPSELFENERGY_H
