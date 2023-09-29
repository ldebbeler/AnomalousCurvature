#include<iostream>
#include "scaling.h"
#include "bubble.h"
#include "timer.hpp"
#include "datastructs.h"
#include "testBubble.h"
#include "writeImag.h"
#include "constants.h"
#include "selfenergy.h"

int main(int argc, char** argv){

    Timer t;

    std::cout << argv[1] << '\n';

    double at{ std::stod(argv[1]) };

    const std::string fileDirectory{ "data/genericOneLoop/data" };
    const std::string suffix{ std::to_string((int)(at*100)) };
    const std::string fileName{ fileDirectory + suffix + ".h5" };

    bubble bub(at);
    //std::cout << "Bubble calculated\n";

    scalingValues scale;
    bubbleValues bv;

    testBubble foo(bub);
    foo.writeBubble(scale, bv);

    selfEnergyValues sev;

    selfenergy bar(at, bub);

    //std::cout << "Positive frequency scaling coefficient: \t" << bar.freqScalePos() << '\n';
    //std::cout << "Negative frequency scaling coefficient: \t" << bar.freqScaleNeg() << '\n';

    double fp{ bar.freqScalePos() };
    double fm{ bar.freqScaleNeg() };
    sev.m_freqPlus = fp;
    sev.m_freqMinus = fm;
    // why does this lead to a segmentation fault?
    //sev.m_freqPlus = bar.freqScalePos();
    //sev.m_freqMinus = bar.freqScaleNeg();

    std::vector<double> freqs(201);
    std::vector<double> selfE(freqs.size());
    for(std::size_t i=0; i<freqs.size(); i++){
        double omega{ -1.0 + 0.01*i };
        freqs[i] = omega;
        selfE[i] = bar.integralLog(omega, 0.0, 0.0).imag();
    }
    sev.m_freqs = freqs;
    sev.m_SE = selfE;

    std::vector<double> rads(201);
    std::vector<double> radPos(rads.size());
    std::vector<double> radNeg(rads.size());
    for(std::size_t i=0; i<rads.size(); i++){
        double x{ -5.0 + 0.05*i };
        rads[i] = x;
        radPos[i] = bar.integralPosx(x);
        radNeg[i] = bar.integralNegx(x);
    }
    sev.m_rads = rads;
    sev.m_radScalingPos = radPos;
    sev.m_radScalingNeg = radNeg;
    
    std::vector<double> tangs(101);
    std::vector<double> tangPos(tangs.size());
    std::vector<double> tangNeg(tangs.size());
    for(std::size_t i=0; i<101; i++){
        double t{ 0.05*i };
        tangs[i] = t;
        tangPos[i] = bar.integralPosTang(t);
        tangNeg[i] = bar.integralNegTang(t);
    }
    sev.m_tangs = tangs;
    sev.m_tangScalingPos = tangPos;
    sev.m_tangScalingNeg = tangNeg;

    H5::H5File file(fileName, H5F_ACC_TRUNC );
    writeh5 writefile(scale, bv, sev);
    writefile.writeMainResults(file);

    std::cout << "Results written to: \n" << fileName << '\n';

    //double tseconds{ t.elapsed() };
    
    //std::cout << "Time taken: " << tseconds << " seconds\n";
 
 

    return 0;
}
