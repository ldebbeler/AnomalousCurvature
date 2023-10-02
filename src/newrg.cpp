#include<iostream>
#include<iomanip>
#include "scaleEqns.h"
#include "datastructs.h"
#include "writeNew.h"
#include "bubbleScale.h"
#include "newSe.h"
#include "radialScale.h"
#include "seRad.h"
#include "seTang.h"
#include "timer.hpp"

int main(int argc, char** argv){
    Timer t;

    //for(int i = 0; i<41; i++){
    int ID{ std::stoi(argv[1]) };
    //int ID = (int)( 400 - 5*i );
    double at{ ID*0.01 };
    //double at{ 2.5 };

    std::string fileID{ std::to_string(ID) };
    std::string name{ prefix + fileID + suffix };
    
    //std::cout << std::setprecision(10) << std::fixed;
/*
    seTang set(at);
    std::cout << "DPlus: \t" << set.m_scalePos.getQuartic() << '\n';
    std::cout << "DMinus: \t" << set.m_scaleNeg.getQuartic() << '\n';

    for(int i=0; i<201; i++){
        double x{ 0.05*i };
        std::cout << '{' << x << ", " << set.m_scalePos.evaluate(x) << "}, \t";
    }
*/

    //newSe foo(at);

    radialScalePos foo(at);
    radialScaleNeg fooNeg(at);

    seRad ser(foo, fooNeg);


    std::cout << "Objects created.\n";

    std::vector<double> x(201);
    std::vector<double> imag(201);
    std::vector<double> real(201);
    std::vector<double> imagDeriv(201);
    std::vector<double> arPos(201);
    std::vector<double> arNeg(201);
    std::vector<double> radInterPos(201);
    std::vector<double> radInterNeg(201);
    for(std::size_t i=0; i<x.size(); i++){
        double y{ -5.0 + i*0.05 };
        x[i] = y;
        imag[i] = foo.bubbleScale(y).imag();
        //imag[i] = foo.exactBubbleQuartic(y).imag();
        real[i] = foo.bubbleScale(y).real();
        //real[i] = foo.exactBubbleQuartic(y).real();
        imagDeriv[i] = foo.imagBubbleDeriv(y);
        arPos[i] = foo.arPos(y);
        arNeg[i] = foo.arNeg(y);
        radInterPos[i] = foo.evaluate(y);
        radInterNeg[i] = fooNeg.evaluate(y);
    }
    std::cout << "Radial Scale Function.\n";

    scalingValues bar;
    bar.m_x = x;
    bar.m_imag = imag;
    bar.m_real = real;
    bar.m_imagDeriv = imagDeriv;
    bar.m_arPos = arPos;
    bar.m_arNeg = arNeg;
    bar.m_apos = foo.Apos();
    bar.m_aneg = foo.Aneg();
    bar.m_radInterPos = radInterPos;
    bar.m_radInterNeg = radInterNeg;
    bar.m_T = at;

    int N{ 50 };
    double a{ 1e-4 };
    double c{ 1.0 };
    double base{ std::exp(std::log(c/a)/(N-1)) };

    std::vector<double> kr(50);
    std::vector<double> sePos(50);
    std::vector<double> seNeg(50);
    for(std::size_t i=0; i<kr.size(); i++){
        double k{ a*std::pow(base,i) };
        kr[i] = k;
        sePos[i] = ser.radValue(k);
        seNeg[i] = ser.radValue(-k);
    }

    std::cout << "Radial Self-Energy.\n";

    std::vector<double> ktilde(501);
    std::vector<double> tangPos(501);
    std::vector<double> tangNeg(501);
    for(std::size_t i=0; i<ktilde.size(); i++){
        double k{ 0.01*i };
        ktilde[i] = k;
        tangPos[i] = foo.atPos(k);
        tangNeg[i] = foo.atNeg(k);
    }

    bar.m_ktilde = ktilde;
    bar.m_atPos = tangPos;
    bar.m_atNeg = tangNeg;

    seTang set(at);

    std::vector<double> interk(501);
    std::vector<double> subtractPos(501);
    std::vector<double> subtractNeg(501);
    std::vector<double> omegas(501);
    std::vector<double> freqPos(501);
    std::vector<double> freqNeg(501);
    for(std::size_t i=0; i<interk.size(); i++){
        double k{ 0.01*i };
        interk[i] = k;
        subtractPos[i] = set.m_scalePos.evaluate(k);
        subtractNeg[i] = set.m_scaleNeg.evaluate(k);
        double omeg{ 0.01*i };
        omegas[i] = omeg;
        freqPos[i] = set.m_scalePos.freqFunction(omeg);
        freqNeg[i] = set.m_scaleNeg.freqFunction(omeg);
    }


    bar.m_interk = interk;
    bar.m_subtractPos = subtractPos;
    bar.m_subtractNeg = subtractNeg;
    bar.m_tangFreqs = omegas;
    bar.m_tangFreqPos = freqPos;
    bar.m_tangFreqNeg = freqNeg;
    bar.m_DPos = set.m_scalePos.getQuartic();
    bar.m_DNeg = set.m_scaleNeg.getQuartic();

    std::cout << "Tangential ScaleFunction.\n";

    std::vector<double> kt(50);
    std::vector<double> seTang(50);
    a = std::pow(10,-2);
    c = 0.1;
    base = std::exp(std::log(c/a)/(N-1));
    for(std::size_t i=0; i<kt.size(); i++){
        double k{ a*std::pow(base,i) };
        kt[i] = k;
        seTang[i] = set.tangValue(k);
    }
    std::cout << "Tangential SelfEnergy.\n";

    bar.m_kr = kr;
    bar.m_kt = kt;
    bar.m_seRadPos = sePos;
    bar.m_seRadNeg = seNeg;
    bar.m_seTang = seTang;

    writeNew write(bar);
    

    H5::H5File file(name, H5F_ACC_TRUNC );
    writeNew writefile(bar);
    writefile.writeMainResults(file);
    
    std::cout << "File created: \t" << name << '\n';
    //}

    double tseconds{ t.elapsed() };
    int minutes{ (int)tseconds/60 };
    int seconds{ (int) tseconds-60*minutes };

    std::cout << "Time taken: " << minutes << " Minutes and " << seconds << " seconds\n";
    return 0;
}
