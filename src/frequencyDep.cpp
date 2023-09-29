#include <iostream>
#include "grid.h"
#include "bubble.h"
#include "interpolate1d.h"
#include "datastructs.h"
#include "selfenergy.h"
#include "testBubble.h"
#include "writeImag.h"
#include "timer.hpp"
#include "constants.h"

int main() {

    Timer t;

    bubble bub;

    scalingValues scale;
    bubbleValues bv;

    testBubble foo(bub);
    foo.writeBubble(scale, bv);

    selfenergy SE(bub); 

    std::cout << SE.integralPosTang(0.1) << '\n';
    std::cout << SE.integralPosx(0.1) << '\n';
    std::cout << SE.integralNegTang(0.1) << '\n';
    std::cout << SE.integralNegx(0.1) << '\n';

    std::vector<double> x(201);
    std::vector<double> funcPos(x.size());
    std::vector<double> funcNeg(x.size());
    
    for(std::size_t i=0; i<x.size(); i++){
        x[i] = -5.0+0.05*i;
        //x[i] = 0.1*i;
    }

    std::cout << x[0] << '\n';
    std::cout << x[x.size()-1] << '\n';
    

    for(std::size_t i=0; i<x.size(); i++){
        funcPos[i] = SE.integralPosx(x[i]);
        funcNeg[i] = SE.integralNegx(x[i]);
    }

    H5::H5File file(filenameFx, H5F_ACC_TRUNC );
    writeh5 writeFile;
    writeFile.writeFunc(file, x, funcPos, funcNeg);


    std::cout << "Results written to: \n " << filenameFx << '\n';

    /*
    std::vector<double> freqs(201);
    std::vector<std::vector<double>> SEValues;
    std::vector<double> kt(21);
    double kr{ 0.0 };
    
    for(std::size_t i=0; i<freqs.size(); i++){
        freqs[i] = -10.0+i*0.1;
    }

    //std::vector<double> krad(3);
    for(std::size_t i=0; i<kt.size(); i++){
        kt[i] = 0.1*i;
    }

    std::cout << freqs[0] << '\n';
    std::cout << freqs[freqs.size()-1] << '\n';
    std::cout << kt[0] << '\n';
    std::cout << kt[kt.size()-1] << '\n';


    for(std::size_t i=0; i<freqs.size(); i++){
        std::vector<double> SEinner(kt.size());
        for(std::size_t j=0; j<kt.size(); j++){
            SEinner[j] = SE.integralLog(freqs[i], kr, kt[j]).imag();
        }
        SEValues.push_back(SEinner);
    }

    std::cout << SEValues.size() << '\n';
    std::cout << SEValues[1].size() << '\n';
    
    seValuesImag sev;
    sev.m_freqs = freqs;
    sev.m_rads = kt;
    sev.m_SE = SEValues;
    sev.m_kr = kr;

    H5::H5File file(filenameImag, H5F_ACC_TRUNC );
    writeh5 writefile(scale, bv, sev);
    writefile.writeMainResults(file);


    std::cout << "Results written to: \n " << filenameImag << '\n';
    */

    double tseconds{ t.elapsed() };
    int minutes = (int) tseconds/60.0;
    int seconds = (int) tseconds - 60.0*minutes;
    
    std::cout << "Time taken: " << minutes << "minutes and " << seconds << " seconds\n";
 
    return 0;
}
