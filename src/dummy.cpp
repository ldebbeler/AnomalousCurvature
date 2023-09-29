#include<iostream>
#include "scaling.h"
#include "bubble.h"
#include "timer.hpp"
#include "datastructs.h"
#include "testBubble.h"
#include "writeImag.h"
#include "constants.h"
#include "indepselfenergy.h"

int main(int argc, char** argv){

    Timer t;

    std::cout << argv[1] << '\n';

    double at{ std::stod(argv[1]) };

    indepselfenergy bar(at);

    const std::string fileDirectory{ "data/genericOneLoop/subtraction" };
    const std::string suffix{ std::to_string((int)(at*100)) };
    const std::string fileName{ fileDirectory + suffix + ".h5" };

    scalingValues scale;
    scale.m_x = bar.getxVector();
    scale.m_real = bar.getReals();
    scale.m_imag = bar.getImags();
    //scale.m_firstReal = bar.getFirstDReal();
    //scale.m_firstImag = bar.getFirstDImag();
    //scale.m_secondReal = bar.getSecondDReal();
    //scale.m_secondImag = bar.getSecondDImag();
    scale.m_T = at;


    selfEnergyValues sev;

    std::cout << "Positive frequency scaling coefficient: \t" << bar.freqScalePos() << '\n';
    std::cout << "Negative frequency scaling coefficient: \t" << bar.freqScaleNeg() << '\n';

    double fp{ bar.freqScalePos() };
    double fm{ bar.freqScaleNeg() };
    sev.m_freqPlus = fp;
    sev.m_freqMinus = fm;

    // frequency dependence of self energy
    std::vector<double> freqs(201);
    std::vector<double> selfE(freqs.size());
    for(std::size_t i=0; i<freqs.size(); i++){
        double omega{ -1.0 + 0.01*i };
        freqs[i] = omega;
        //selfE[i] = bar.integralLog(omega, 0.0, 0.0).imag();
    }
    sev.m_freqs = freqs;
    sev.m_SE = selfE;

    // radial momentum dependence of self energy
    std::vector<double> rads(201);
    std::vector<double> radPos(rads.size());
    std::vector<double> radNeg(rads.size());
    for(std::size_t i=0; i<rads.size(); i++){
        double x{ -5.0 + 0.05*i };
        rads[i] = x;
        //radPos[i] = bar.integralPosx(x);
        //radNeg[i] = bar.integralNegx(x);
    }
    sev.m_rads = rads;
    sev.m_radScalingPos = radPos;
    sev.m_radScalingNeg = radNeg;
    
    // tangential momentum dependence of self energy
    std::vector<double> tangs(301);
    std::vector<double> tangPos(tangs.size());
    std::vector<double> tangNeg(tangs.size());
    for(std::size_t i=0; i<tangs.size(); i++){
        double t{ 0.001*(i+1) };
        tangs[i] = t;
        tangPos[i] = bar.integralPosTang(t);
        tangNeg[i] = bar.integralNegTang(t);
    }
    sev.m_tangs = tangs;
    sev.m_tangScalingPos = tangPos;
    sev.m_tangScalingNeg = tangNeg;

    /*
    std::vector<double> dt(101);
    std::vector<double> difftangPos(dt.size());
    std::vector<double> difftangNeg(dt.size());
    std::vector<double> derivPos(dt.size());
    std::vector<double> derivNeg(dt.size());

    for(std::size_t i=0; i<101; i++){
        double t{ 0.005*i };
        dt[i] = t;
        difftangPos[i] = bar.twoDPos(t);
        difftangNeg[i] = bar.twoDNeg(t);
        //derivPos[i] = bar.derivPosValue(t);
        //derivNeg[i] = bar.derivNegValue(t);
    }
    sev.m_dkt= dt;
    sev.m_dTangPos = difftangPos;
    sev.m_dTangNeg = difftangNeg;
    sev.m_derivKt = dt;
    sev.m_derivTangPos = derivPos;
    sev.m_derivTangNeg = derivNeg;
*/

    H5::H5File file(fileName, H5F_ACC_TRUNC );
    writeh5 writefile(scale, sev);
    writefile.writeMainResults(file);

    std::cout << "Results written to: \n" << fileName << '\n';

    double tseconds{ t.elapsed() };
    int minutes = (int) tseconds/60.0;
    int seconds = (int) tseconds - 60.0*minutes;
    std::cout << "Time taken: " << minutes << " minutes and " << seconds << " seconds\n";
 
    return 0;
}
