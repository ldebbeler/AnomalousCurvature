#include "selfenergy.h"
#include "constants.h"
#include "d2.hpp"
#include "d3.hpp"
#include "integrate2d.hpp"
#include "integrate2d_dcuhre.hpp"
#include "integrate3d.hpp"

double selfenergy::dispersion(double qr, double qt){
    return vF*qr + b*std::pow(std::abs(qt),4);
}

double selfenergy::dispersion(double kr, double qr, double kt, double qt){
    return -vF*(kr-qr) + b*std::pow(std::abs(kt-qt),4);
}

selfenergy::selfenergy(bubble foo): m_bubble(foo) {}

std::complex<double> selfenergy::integrandPos(double omega, double kr, double kt, double qr, double qt){
    double om{ std::abs(omega) };
    return 1.0/m_bubble.evaluateCombo(vF*qr-om,qr-b/vF*std::pow(qt-kt,4)+kr,qt);
}

std::complex<double> selfenergy::integrandPosLog(double omega, double kr, double kt, double qr, double Q){
    return std::exp(Q)*(integrandPos(omega,kr,kt,qr,std::exp(Q)) + integrandPos(omega,kr,kt,qr,-std::exp(Q)));
}

std::complex<double> selfenergy::integrandNeg(double omega, double kr, double kt, double qr, double qt){
    double om{ std::abs(omega) };
    return -1.0/m_bubble.evaluateCombo(vF*qr+om,qr-b/vF*std::pow(qt-kt,4)+kr,qt);
}

std::complex<double> selfenergy::integrandNegLog(double omega, double kr, double kt, double qr, double Q){
    return std::exp(Q)*(integrandNeg(omega,kr,kt,qr,std::exp(Q)) + integrandNeg(omega,kr,kt,qr,-std::exp(Q)));
}

// no factor 2 here because integral w.r.t. qt not symmetric
std::complex<double> selfenergy::integralLog(double omega, double kr, double kt){
    if(omega>0.0){
        auto integ = [&] (double2 q) -> std::complex<double>{
            double qr = q.x;
            double Q = q.y;
            return integrandPosLog(omega, kr, kt, qr, Q);
        };
        //return integrate2Dmeasure2pi(0.0, omega/vF, std::log(IR), std::log(UV)).integrate(integ, prec2d, steps2d);
        return averageBZ(0., omega/vF, std::log(IR), std::log(UV)).integrate_dcuhre(integ, prec2d);
    }
    else{
        auto integ = [&] (double2 q) -> std::complex<double>{
            double qr = q.x;
            double Q = q.y;
            return integrandNegLog(omega, kr, kt, qr, Q);
        };
        return integrate2Dmeasure2pi(omega/vF, 0.0, std::log(IR), std::log(UV)).integrate(integ, prec2d, steps2d);
    }
}

std::complex<double> selfenergy::integrandPosx(double qr, double qt, double x){
    return -1.0/m_bubble.evaluateCombo(vF*qr-1.0,qr+x-b/vF*std::pow(qt,4),qt);
}

//factor of 2 because qt integral symmetric
std::complex<double> selfenergy::integrandPosLogx(double qr, double Qt, double x){
    return 2.0*std::exp(Qt)*integrandPosx(qr, std::exp(Qt), x);
}

double selfenergy::integralPosx(double x){
    auto integ = [&] (double2 q) -> std::complex<double>{
        double qr = q.x;
        double Qt = q.y;
        return integrandPosLogx(qr, Qt, x);
    };
    return integrate2Dmeasure2pi(0.0, 1.0/vF, std::log(IR), std::log(UV)).integrate(integ, prec2d, steps2d).imag();
}

std::complex<double> selfenergy::integrandNegx(double qr, double qt, double x){
    return 1.0/m_bubble.evaluateCombo(vF*qr+1.0,qr+x-b/vF*std::pow(qt,4),qt);
}

std::complex<double> selfenergy::integrandNegLogx(double qr, double Qt, double x){
    return 2.0*std::exp(Qt)*integrandNegx(qr, std::exp(Qt), x);
}

double selfenergy::integralNegx(double x){
    auto integ = [&] (double2 q) -> std::complex<double>{
        double qr = q.x;
        double Qt = q.y;
        return integrandNegLogx(qr, Qt, x);
    };
    return integrate2Dmeasure2pi(-1.0/vF, 0.0, std::log(IR), std::log(UV)).integrate(integ, prec2d, steps2d).imag();
    //return averageBZ(-1.0/vF, 0.0, std::log(IR), std::log(UV)).integrate_dcuhre(integ, prec2d);
}



std::complex<double> selfenergy::integrandPosTang(double qr, double qt, double x){
    return -1.0/m_bubble.evaluateCombo(vF*qr-1.0,qr-b/vF*std::pow(qt-x,4),qt);
}

std::complex<double> selfenergy::integrandPosLogTang(double qr, double Qt, double x){
    return std::exp(Qt)*(integrandPosTang(qr, std::exp(Qt), x)+integrandPosTang(qr,-std::exp(Qt), x));
}

double selfenergy::integralPosTang(double x){
    auto integ = [&] (double2 q) -> std::complex<double>{
        double qr = q.x;
        double Qt = q.y;
        return integrandPosLogTang(qr, Qt, x) - integrandPosLogTang(qr, Qt, 0.0);
    };
    //return integrate2Dmeasure2pi(0.0, 1.0/vF, std::log(IR), std::log(UV)).integrate(integ, prec2d, steps2d).imag();
    return averageBZ(0.0, 1.0/vF, std::log(IR), std::log(UV)).integrate_dcuhre(integ, prec2d).imag();
}

std::complex<double> selfenergy::integrandNegTang(double qr, double qt, double x){
    return 1.0/m_bubble.evaluateCombo(vF*qr+1.0,qr-b/vF*std::pow(qt-x,4),qt);
}

std::complex<double> selfenergy::integrandNegLogTang(double qr, double Qt, double x){
    return std::exp(Qt)*(integrandNegTang(qr, std::exp(Qt), x) + integrandNegTang(qr,-std::exp(Qt), x));
}

double selfenergy::integralNegTang(double x){
    auto integ = [&] (double2 q) -> std::complex<double>{
        double qr = q.x;
        double Qt = q.y;
        return integrandNegLogTang(qr, Qt, x) - integrandNegLogTang(qr, Qt, 0.0);
    };
    //return integrate2Dmeasure2pi(-1.0/vF, 0.0, std::log(IR), std::log(UV)).integrate(integ, prec2d, steps2d).imag();
    return averageBZ(-1.0/vF, 0.0, std::log(IR), std::log(UV)).integrate_dcuhre(integ, prec2d).imag();
}




// Real part *************************************************************************************
// kr dependence
std::complex<double> selfenergy::fluctuationPropagator(double nu, double qr, double qt){
    return -1.0/(m_bubble.evaluateCombo(nu, qr, qt));
}

std::complex<double> selfenergy::BoseIntegrand1(double nu, double qr, double qt, double kr){
    return 4.0*fluctuationPropagator(nu, qr+kr+nu/vF -b/vF*std::pow(qt,4), qt).imag()/(-vF*qr);
}

std::complex<double> selfenergy::BoseIntegrand2(double nu, double qr, double qt, double kr){
    return 4.0*fluctuationPropagator(nu,-qr+kr+nu/vF -b/vF*std::pow(qt,4), qt).imag()/( vF*qr);
}

std::complex<double> selfenergy::FermiIntegrand(double nu, double qr, double qt, double kr){
    return -4.0*M_PI/(UV-IR)*fluctuationPropagator(vF*qr, qr+kr-b/vF*std::pow(qt,4), qt).real();
}

std::complex<double> selfenergy::RadialIntegrand(double nu, double qr, double qt, double kr){
    return BoseIntegrand1(nu,qr,qt,kr)-BoseIntegrand1(nu,qr,qt,0.0)
          +BoseIntegrand2(nu,qr,qt,kr)-BoseIntegrand2(nu,qr,qt,0.0)
          +FermiIntegrand(nu,qr,qt,kr)-FermiIntegrand(nu,qr,qt,0.0);
}

std::complex<double> selfenergy::RadialLogIntegrand(double Qn, double Qr, double Qt, double kr){
    return std::exp(Qn)*std::exp(Qr)*std::exp(Qt)*(RadialIntegrand(-std::exp(Qn),-std::exp(Qr),-std::exp(Qt), kr));
}

double selfenergy::radialValue(double kr){
    double ir{ std::log(IR) };
    double uv{ std::log(UV) };
    auto integ = [&] (double3 Q) -> std::complex<double>{
        double Qn = Q.x;
        double Qr = Q.y;
        double Qt = Q.z;
        return RadialLogIntegrand(Qn, Qr, Qt, kr);
    };
    double radValue = integrate3Dmeasure2pi(ir,uv,ir,uv,ir,uv).integrate(integ, prec3d, steps3d).real();
    return radValue;
}

/*
std::complex<double> selfenergy::greenFunction(double omega, double kr, double kt){
    return 1.0/(omega + vF*kr - b*std::pow(kt,4));
}

// factor 2 here, to account for 1/2pi while in reality we have only 1/pi
std::complex<double> selfenergy::realIntegrand1(double omega, double kr, double kt, double nu, double qr, double qt){
    return 2.0*fluctuationPropagator(nu,qr-kr,qt-kt).imag()*greenFunction(omega+nu,qr,qt)
        -  2.0*fluctuationPropagator(nu,qr,qt).imag()*greenFunction(nu,qr,qt);
}

std::complex<double> selfenergy::realIntegrand2(double omega, double kr, double kt, double qr, double qt){
    return fluctuationPropagator(-vF*qr-omega, qr+b/vF*std::pow(qt,4)-kr,qt-kt).real()
          -fluctuationPropagator(-vF*qr,qr+b/vF*std::pow(qt,4),qt).real();
}

std::complex<double> selfenergy::realLogIntegrand1(double omega, double kr, double kt, double Qn, double Qr, double Qt){
    return std::exp(Qn)*std::exp(Qr)*std::exp(Qt)*(realIntegrand1(omega, kr, kt, std::exp(Qn), std::exp(Qr), std::exp(Qt))
                                                  +realIntegrand1(omega, kr, kt, std::exp(Qn),-std::exp(Qr), std::exp(Qt))
                                                  +realIntegrand1(omega, kr, kt, std::exp(Qn), std::exp(Qr),-std::exp(Qt))
                                                  +realIntegrand1(omega, kr, kt, std::exp(Qn),-std::exp(Qr),-std::exp(Qt)));
}

std::complex<double> selfenergy::realLogIntegrand2(double omega, double kr, double kt, double Qr, double Qt){
    return std::exp(Qr)*std::exp(Qt)*(realIntegrand2(omega, kr, kt, std::exp(Qr), std::exp(Qt))
                                     +realIntegrand2(omega, kr, kt, std::exp(Qr),-std::exp(Qt)));
}

std::complex<double> selfenergy::realLogIntegral1(double omega, double kr, double kt){
    double ir{ std::log(IR) };
    double uv{ std::log(UV) };
    double irt{ std::log(IRt) };
    double uvt{ std::log(UVt) };
    auto integ = [&] (double3 Q) -> std::complex<double>{
        double Qn = Q.x;
        double Qr = Q.y;
        double Qt = Q.z;
        return realLogIntegrand1(omega, kr, kt, Qn, Qr, Qt) + realLogIntegrand2(omega, kr, kt, Qr, Qt)*2.0*M_PI/(uv-ir);
    };
    double realValue1 = integrate3Dmeasure2pi(ir,uv,ir,uv,irt,uvt).integrate(integ,prec3d,steps3d).real();
    return realValue1;
}

std::complex<double> selfenergy::realLogIntegral2(double omega, double kr, double kt){
    auto integ = [&] (double2 Q) -> std::complex<double>{
        double Qr = Q.x;
        double Qt = Q.y;
        return realLogIntegrand2(omega, kr, kt, Qr, Qt);
    };
    double ir{ std::log(IR)  };
    double uv{ std::log(UV)  };
    double realValue2 = integrate2Dmeasure2pi(ir,uv,ir,uv).integrate(integ,prec2d,steps2d).real();
    return realValue2;
}
*/
/*
double selfenergy::realValue(double omega, double kr, double kt){
    return realLogIntegral1(omega, kr, kt).real() + realLogIntegral2(omega, kr, kt).real();
}
*/
