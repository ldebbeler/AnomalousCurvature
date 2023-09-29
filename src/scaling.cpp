#include "scaling.h"
//#include "integrategsl.hpp"
#include "constants.h"
#include "d2.hpp"
#include "integrate2d.hpp"

scaling::scaling(double T): m_T(T) {}

std::complex<double> scaling::termx(double k0, double kt, double x){
    return 1.0/(2.0*I*k0-std::pow(std::abs(kt-0.5),m_T) - std::pow(std::abs(kt+0.5),m_T) + x);
}

std::complex<double> scaling::termZero(double k0, double kt){
    return 1.0/(2.0*I*k0 - 2*std::pow(std::abs(kt),m_T));
}

std::complex<double> scaling::integrand(double k0, double kt, double x){
    return -4.0*I*(termx(k0,kt,x)-termZero(k0,kt));
}

std::complex<double> scaling::logIntegrand(double K, double Kt, double x){
    return std::exp(K)*std::exp(Kt)*(integrand(std::exp(K), std::exp(Kt), x) + integrand(std::exp(K),-std::exp(Kt), x));
}

std::complex<double> scaling::value(double x){
    auto integ = [&] (double2 k) -> std::complex<double>{
        double K = k.x;
        double Kt = k.y;
        return logIntegrand(K, Kt, x);
    };
    return integrate2Dmeasure2pi(std::log(IR), std::log(UV), std::log(IR), std::log(UV)).integrate(integ, prec2d, steps2d);
}

std::complex<double> scaling::firstDIntegrand(double k0, double kt, double x){
    return 4.0*I/(std::pow(2.0*I*k0-std::pow(std::abs(kt-0.5),m_T)-std::pow(std::abs(kt+0.5),m_T) + x,2.0));
}

std::complex<double> scaling::firstDlogIntegrand(double K0, double Kt, double x){
    return std::exp(K0)*std::exp(Kt)*(firstDIntegrand(std::exp(K0), std::exp(Kt), x) + firstDIntegrand(std::exp(K0),-std::exp(Kt), x));
}

std::complex<double> scaling::firstDvalue(double x){
    auto integ = [&] (double2 k) -> std::complex<double>{
        double K0 = k.x;
        double Kt = k.y;
        return firstDlogIntegrand(K0, Kt, x);
    };
    return integrate2Dmeasure2pi(std::log(IR), std::log(UV), std::log(IR), std::log(UV)).integrate(integ, prec2d, steps2d);
}

std::complex<double> scaling::secondDIntegrand(double k0, double kt, double x){
    return -8.0*I/(std::pow(2.0*I*k0-std::pow(std::abs(kt-0.5),m_T)-std::pow(std::abs(kt+0.5),m_T) + x,3.0));
}

std::complex<double> scaling::secondDlogIntegrand(double K0, double Kt, double x){
    return std::exp(K0)*std::exp(Kt)*(secondDIntegrand(std::exp(K0), std::exp(Kt), x) + secondDIntegrand(std::exp(K0),-std::exp(Kt), x));
}

std::complex<double> scaling::secondDvalue(double x){
    auto integ = [&] (double2 k) -> std::complex<double>{
        double K0 = k.x;
        double Kt = k.y;
        return secondDlogIntegrand(K0, Kt, x);
    };
    return integrate2Dmeasure2pi(std::log(IR), std::log(UV), std::log(IR), std::log(UV)).integrate(integ, prec2d, steps2d);
}



std::complex<double> scaling::coefficientZeroIntegrand(double k0, double kt, double sx){
    return -4.0*I*(1.0/(2.0*I*k0 - 2.0*std::pow(std::abs(kt),m_T) + sx) - 1.0/(2.0*I*k0 - 2.0*std::pow(std::abs(kt),m_T)));
}

std::complex<double> scaling::coefficientZero(double sx){
    auto integ = [&] (double2 k) -> std::complex<double>{
        double K = k.x;
        double Kt = k.y;
        return std::exp(K)*std::exp(Kt)*(coefficientZeroIntegrand(std::exp(K), std::exp(Kt), sx)
                                        +coefficientZeroIntegrand(std::exp(K),-std::exp(Kt), sx));
    };
    return integrate2Dmeasure2pi(std::log(IR), std::log(UV), std::log(IR), std::log(UV)).integrate(integ, prec2d, steps2d);
}

std::complex<double> scaling::coefficientTwoIntegrand(double k0, double kt, double sx){
    return (-I*m_T*(m_T-1.0)*std::pow(std::abs(kt),m_T-2.0))/(std::pow(2.0*I*k0-2*std::pow(std::abs(kt),m_T)+sx,2.0));
}

std::complex<double> scaling::coefficientTwo(double sx){
    auto integ = [&] (double2 k) -> std::complex<double>{
        double K = k.x;
        double Kt = k.y;
        return std::exp(K)*std::exp(Kt)*(coefficientTwoIntegrand(std::exp(K), std::exp(Kt), sx)
                                        +coefficientTwoIntegrand(std::exp(K),-std::exp(Kt), sx));
    };
    return integrate2Dmeasure2pi(std::log(IR), std::log(UV), std::log(IR), std::log(UV)).integrate(integ, prec2d, steps2d);
}
