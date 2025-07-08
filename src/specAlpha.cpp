#include<iostream>
#include<iomanip>
#include "scaleEqns.h"
#include "datastructs.h"
#include "writeFile.h"
#include "bubbleScale.h"
#include "selfenergy.h"
#include "radialScale.h"
#include "seRad.h"
#include "tangentialScale.h"
#include "seTang.h"
#include "timer.hpp"

int main(int argc, char** argv){
    Timer t;

    //pass alpha*100 as command line input
    int ID{ std::stoi(argv[1]) };
    double at{ ID*0.01 };

    std::string fileID{ std::to_string(ID) };
    std::string name{ prefix + fileID + suffix };
    
    //create self energy object
    selfenergy se(at);
    scaleEqns bubbleI(at);
    // scaling form of bubble
    bubbleScaleReal realBubble(at);
    bubbleScaleImag imagBubble(at);
    // calculate scaling function
    std::vector<double> x_bscale(1001);
    std::vector<double> real_bscale(1001);
    std::vector<double> imag_bscale(1001);
    for(std::size_t i=0; i<x_bscale.size(); i++){
        double x{ -5.001 + i*0.01 };
        x_bscale[i] = x;
        imag_bscale[i] = bubbleI.imagValue(x);
        real_bscale[i] = bubbleI.realValue(x);
    }

    // frequency coefficient of self energy
    double Apos{ se.Apos() };
    double Aneg{ se.Aneg() };


    // radial momentum scaling function of self energy
    std::vector<double> krt(101);
    std::vector<double> arPos(101);
    std::vector<double> arNeg(101);
    for(std::size_t i=0; i<krt.size(); i++){
        double arg{ -5.0 + 0.1*i };
        krt[i] = arg;
        arPos[i] = se.arPos(arg);
        arNeg[i] = se.arNeg(arg);
    }

    // radial momentum dependence via Kramers-Kronig
    seRad radial(at);
    std::vector<double> kr(100);
    std::vector<double> seRad(100);
    for(std::size_t i=0; i<kr.size(); i++){
        double arg{ 0.01*(i+1) };
        kr[i] = arg;
        seRad[i] = radial.radValue(arg);
    }

    // tangential momentum scaling function of self energy
    std::vector<double> ktt(101);
    std::vector<double> atPos(101);
    std::vector<double> atNeg(101);
    for(std::size_t i=0; i<ktt.size(); i++){
        double arg{ 0.05*i };
        ktt[i] = arg;
        atPos[i] = se.atPos(arg);
        atNeg[i] = se.atNeg(arg);
    }

    // tangential momentum dependence via Kramers-Kronig
    seTang tangential(at);
    std::vector<double> kt(100);
    std::vector<double> seTang(100);
    for(std::size_t i=0; i<kt.size(); i++){
        double arg{ 0.01*(i+1) };
        kt[i] = arg;
        seTang[i] = tangential.tangValue(arg);
    }

    double deltab = tangential.deltab();

    // store results in datastructs
    scalingValues results;

    results.m_T = at;

    results.m_x = x_bscale;
    results.m_real = real_bscale;
    results.m_imag = imag_bscale;

    results.m_apos = Apos;
    results.m_aneg = Aneg;

    results.m_krt = krt;
    results.m_arPos = arPos;
    results.m_arNeg = arNeg;

    results.m_kr = kr;
    results.m_seRad = seRad;

    results.m_ktt = ktt;
    results.m_atPos = atPos;
    results.m_atNeg = atNeg;

    results.m_kt = kt;
    results.m_seTang = seTang;

    results.m_deltab = deltab;
    results.m_Cpos = tangential.m_Cpos;
    results.m_Cneg = tangential.m_Cneg;
    results.m_Dpos = tangential.m_Dpos;
    results.m_Dneg = tangential.m_Dneg;

    H5::H5File file(name, H5F_ACC_TRUNC );
    writeNew writefile(results);
    writefile.writeMainResults(file);
    
    std::cout << "File created: \t" << name << '\n';

    double tseconds{ t.elapsed() };
    int minutes{ (int)tseconds/60 };
    int seconds{ (int) tseconds-60*minutes };

    std::cout << "Time taken: " << minutes << " Minutes and " << seconds << " seconds\n";
    return 0;
}

