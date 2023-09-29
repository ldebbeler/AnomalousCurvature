#include <iostream>
#include "scaling.h"
#include "grid.h"
#include "bubble.h"
#include "interpolate1d.h"
#include "datastructs.h"
#include "selfenergy.h"
#include "testBubble.h"
#include "writeReal.h"
#include "timer.hpp"
#include "constants.h"
#include "omp.h"

//ask pietro/deme whether this is a good idea
//requires func like objects, here with one dependency.
//could be generalized to more arguments perhaps via d2 and d3
template <typename T>
double logFunc(T func, double K){
    return std::exp(K)*(func(std::exp(K))+func(-std::exp(K)));
}

int main() {

    // save bubble values in these structs
    scalingValues scale;
    bubbleValues bv; 
    bubble bub;
    std::cout << "Bubble calculated! \n";
    testBubble foo(bub);
    foo.writeBubble(scale, bv);

    Timer t;
    selfenergy SE(bub);

    H5::H5File file(filenameRealParallel, H5F_ACC_TRUNC );


    std::vector<double> kr(M);
    std::vector<double> seRadialPos(M);
    std::vector<double> seRadialNeg(M);
    double x{ std::exp(std::log(maxk/mink)/(M-1)) };
    //parallelize this loop
    #pragma omp parallel for schedule(dynamic)
    for(std::size_t i=0; i<kr.size(); i++){
        double val{ mink*std::pow(x,i) };
        //std::cout << "Num threads: " << omp_get_num_threads() << std::endl;
        kr[i] = val;
        seRadialPos[i] = SE.radialValue(val);
    }

    seValuesReal sev;
    sev.m_val = kr;
    sev.m_SEPos = seRadialPos;
    sev.m_SENeg = seRadialNeg;

    writeh5 writefile(scale, bv, sev);
    writefile.writeMainResults(file);

    std::cout << "Results written to: \n " << filenameReal << '\n';

    double tseconds{ t.elapsed() };
    int minutes = (int) tseconds/60.0;
    int seconds = (int) tseconds - 60.0*minutes;
    int hours = (int) minutes/60.0;
    minutes = minutes - 60.0*hours;
    
    std::cout << "Time taken: " << hours << " hours, " << minutes << " minutes and " << seconds << " seconds\n";
    std::cout << "Program executed!\n";

    return 0;
}
