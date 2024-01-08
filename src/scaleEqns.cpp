#include "scaleEqns.h"
#include "rootFinder.hpp"
#include "integrategsl.hpp"
#include<iostream>

scaleEqns::scaleEqns(double T): m_T(T) {}


double scaleEqns::rootEqn(double kt, double x){
    return std::pow(std::abs(kt+0.5),m_T) + std::pow(std::abs(kt-0.5),m_T) - x;
}

double scaleEqns::ktzero(double x){
    auto func = [&] (double kt) -> double{ return rootEqn(kt, x); };
    double kt0 = rootFinder().findRoot(func, x_down, x_up, x_down);
    return kt0;
}

double scaleEqns::imagValue(double x){
    if(x>std::pow(2,1-m_T)+ 0.001){
        return -ktzero(x)/M_PI;
    }
    else{
        return 0.0;
    }
}

double scaleEqns::imagAsymptotic(double x){
    if(x<0){ return 0.0; }
    else{
        return -std::pow(std::abs(x),1.0/m_T)/(M_PI*std::pow(2,1.0/m_T));
    }
}

interpolater1d scaleEqns::imagSpline(){
    std::vector<double> x(nodes);
    std::vector<double> imags(nodes);
    double k{ 2*Lambda/(nodes-1) };
    for(int i=0; i<nodes; i++){
        double val{ -Lambda + k*i };
        x[i] = val;
        imags[i] = imagValue(val);
    }
    interpolater1d foo(x,imags);
    return foo;
}

double scaleEqns::realIntegrand(double kt, double x){
    return 1.0/(2.0*M_PI*M_PI)*std::log(std::abs(std::pow(std::abs(kt+0.5),m_T) + std::pow(std::abs(kt-0.5),m_T) - x)/(2*std::pow(std::abs(kt),m_T)));
}

double scaleEqns::realValue(double x){
    // factor 2 because of symmetry w.r.t. kt
    if(x<std::pow(2,1-m_T)+eps){
        auto integ1 = [&] (double K) -> double{ return 2.0*std::exp(K)*realIntegrand(std::exp(K),x); };
        double val1{ integrater().integrate(integ1, std::log(IR), std::log(UV)) };
        return val1;
    }
    else{
        double kzero{ ktzero(x) };
        auto integ1 = [&] (double K) -> double{ return 2.0*std::exp(K)*(realIntegrand(std::exp(K),x) + realIntegrand(kzero-std::exp(K),x)); };
        double val1{ integrater().integrate(integ1, std::log(IR), std::log(0.5*kzero)) };
        auto integ2 = [&] (double K) -> double{ return 2.0*std::exp(K)*(realIntegrand(kzero + std::exp(K),x)); };
        double val2{ integrater().integrate(integ2, std::log(IR), std::log(UV)) };
        return val1 + val2;
    }
}

double scaleEqns::realAsymptotic(double x){
    if(x>0.0){
        return std::pow(std::abs(x),1.0/m_T)/(M_PI*std::pow(2,1.0/m_T))*cos(M_PI/m_T)/sin(M_PI/m_T);
    }
    else{
        return std::pow(std::abs(x),1.0/m_T)/(M_PI*std::pow(2,1.0/m_T))*1.0/sin(M_PI/m_T);
    }
}

interpolater1d scaleEqns::realSpline(){
    std::vector<double> x(nodes);
    std::vector<double> reals(nodes);
    double k{ 2*Lambda/(nodes-1) };
    for(int i=0; i<nodes; i++){
        double val{ -Lambda + k*i };
        x[i] = val;
        reals[i] = realValue(val);
    }
    interpolater1d foo(x,reals);
    return foo;
}


