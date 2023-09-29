#include "grid.h"
#include "constants.h"
#include <cmath>


grid::grid(double init, double fin, std::size_t L): m_x(logSpacing?logCreate(init, fin, L):linCreate(fin,L)) {}

std::vector<double> grid::linCreate(double endpoint, std::size_t L){
    double x{ 2.0*endpoint/(L-1.0) };
    std::vector<double> vec(L);
    for(std::size_t i=0; i<L; i++){
        vec[i] = -endpoint + i*x;
    }
    return vec;
}


std::vector<double> grid::logCreate(double init, double fin, std::size_t L){
    std::vector<double> vec(2*L + 1);
    vec[L]=0.0;
    double x{ std::exp(std::log(fin/init)/(L-1.0)) };
    for (std::size_t i=0; i<L; i++){
        vec[L-1-i] = -init*std::pow(x,i);
        vec[L+1+i] =  init*std::pow(x,i);
    }
    return vec;
}


