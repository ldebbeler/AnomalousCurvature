#include "bubble.h"
#include "grid.h"
#include "constants.h"
#include <iostream>
#include "d2.hpp"
#include "integrate2d.hpp"

bubble::bubble(double T): m_T(T), m_scale(scaling(m_T)), 
                        m_extraPosZero(m_scale.coefficientZero(1.0)), m_extraNegZero(m_scale.coefficientZero(-1.0)),
                        m_extraPosTwo(m_scale.coefficientTwo(1.0)), m_extraNegTwo(m_scale.coefficientTwo(-1.0)),
                        m_x(gridVector()), m_scalingReal(realVector()), m_scalingImag(imagVector()),
                        m_splineReal(interpolater1d(m_x,m_scalingReal)), 
                        m_splineImag(interpolater1d(m_x,m_scalingImag)) {}

std::vector<double> bubble::gridVector(){
    grid vec(minx, extraPolate, N);
    return vec.m_x;
}

std::vector<double> bubble::realVector(){
    std::vector<double> itx(m_x.size());
    for(std::size_t i=0; i<m_x.size(); i++){
        //itx[i] = realLogIntegral(m_x[i]);
        itx[i] = m_scale.value(m_x[i]).real();
    }
    return itx;
}

std::vector<double> bubble::imagVector(){
    std::vector<double> itx(m_x.size());
    for(std::size_t i=0; i<m_x.size(); i++){
        //itx[i] = imagLogIntegral(m_x[i]);
        itx[i] = m_scale.value(m_x[i]).imag();
    }
    return itx;
}
/*
std::complex<double> bubble::quadraticPos(double omega, double qr, double qt){
    double prefactor{ (3*std::pow(2*b,0.25))/(16*M_PI*vF) };
    return prefactor/std::pow(std::abs(omega - vF*qr),0.25)*std::pow(std::abs(qt),2)*(1.0+I)/2.0;
}

std::complex<double> bubble::quadraticNeg(double omega, double qr, double qt){
    double prefactor{ (3*std::pow(2*b,0.25))/(16*M_PI*vF) };
    return prefactor/std::pow(std::abs(omega - vF*qr),0.25)*std::pow(std::abs(qt),2)*1.0/std::sqrt(2.0);
}
*/
// Form for "Plus" Bubble. Other bubble is simply obtained by known relationship: negative freq and complex conjugation
std::complex<double> bubble::evaluatePlus(double omega, double qr, double qt){
    double val{ (omega-vF*qr)/(b*std::pow(std::abs(qt),m_T)) };
    if(val < -extraPolate+0.1){
        return std::abs(qt)/(4*vF)*(std::pow(std::abs(val),1.0/m_T)*m_extraNegZero + std::pow(std::abs(val),-1.0/m_T)*m_extraNegTwo);
    }
    else if(val > extraPolate-0.1){
        return std::abs(qt)/(4*vF)*(std::pow(std::abs(val),1.0/m_T)*m_extraPosZero + std::pow(std::abs(val),-1.0/m_T)*m_extraPosTwo);
    }
    else{
        return std::abs(qt)/(4*vF)*(m_splineReal.evaluate(val)+I*m_splineImag.evaluate(val));
    }
}

std::complex<double> bubble::evaluate(double omega, double qr, double qt){
    return evaluatePlus(omega, qr, qt) + std::conj(evaluatePlus(-omega, qr, qt));
}
/*
std::complex<double> bubble::evaluateStatic(double2 q){
    return evaluateSingle(0.0, q.x, q.y) + std::conj(evaluateSingle(0.0, q.x, q.y));
}
*/
std::vector<double> bubble::getxVector(){
    return m_x;
}

std::vector<double> bubble::getReals(){
    return m_scalingReal;
}

std::vector<double> bubble::getImags(){
    return m_scalingImag;
}

double bubble::getExp(){
    return m_T;
}
