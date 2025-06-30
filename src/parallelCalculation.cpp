#include<iostream>
#include "scaleEqns.h"
#include "datastructs.h"
#include "writeNew.h"
#include "bubbleScale.h"
#include "newSe.h"
#include "radialScale.h"
#include "seRad.h"
#include "seTang.h"
#include "timer.hpp"
#include "omp.h"

int main(){
    Timer t;

    #pragma omp parallel for schedule(dynamic)
    for(int i = 0; i<5; i++){

        int ID = (int)( 275 - 5*i );
        double at{ ID*0.01 };
        //double at{ 2.5 };

        std::string fileID{ std::to_string(ID) };
        std::string name{ prefix + fileID + suffix };
        //std::string name{ "newdata/popel250.h5" };

        //newSe foo(at);

        radialScalePos foo(at);
        radialScaleNeg fooNeg(at);

        seRad ser(foo, fooNeg);

        seTang set(at);

        //std::cout << "Objects created.\n";

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
            real[i] = foo.bubbleScale(y).real();
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
        double b{ 1.0 };
        double base{ std::exp(std::log(b/a)/(N-1)) };

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

        std::vector<double> ktilde(101);
        std::vector<double> tangPos(101);
        std::vector<double> tangNeg(101);
        for(std::size_t i=0; i<ktilde.size(); i++){
            double k{ 0.1*i };
            ktilde[i] = k;
            tangPos[i] = foo.atPos(k);
            tangNeg[i] = foo.atNeg(k);
        }

        bar.m_ktilde = ktilde;
        bar.m_atPos = tangPos;
        bar.m_atNeg = tangNeg;

        std::vector<double> interk(101);
        std::vector<double> subtractPos(101);
        std::vector<double> subtractNeg(101);
        for(std::size_t i=0; i<interk.size(); i++){
            double k{ 0.075*i };
            interk[i] = k;
            subtractPos[i] = set.m_scalePos.evaluate(k);
            subtractNeg[i] = set.m_scaleNeg.evaluate(k);
        }

        bar.m_interk = interk;
        bar.m_subtractPos = subtractPos;
        bar.m_subtractNeg = subtractNeg;

        std::cout << "Tangential ScaleFunction.\n";

        std::vector<double> kt(50);
        std::vector<double> seTang(50);
        b = 0.1;
        base = std::exp(std::log(b/a)/(N-1));
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
    }

    /*
    std::vector<double> alphat(40);
    std::vector<double> apos(40);
    std::vector<double> aneg(40);

    for(std::size_t i=0; i<alphat.size(); i++){
        double at{ 4.0 - 0.05*i };
        newSe foo(at);
        alphat[i] = at;
        apos[i] = foo.Apos();
        aneg[i] = foo.Aneg();

        std::vector<double> x(201);
        std::vector<double> imag(201);
        std::vector<double> real(201);
        for(std::size_t i=0; i<x.size(); i++){
            double y{ -5.0 + i*0.05 };
            x[i] = y;
            imag[i] = foo.bubbleScale(y).imag();
            real[i] = foo.bubbleScale(y).real();
        }

        scalingValues bar;
        bar.m_x = x;
        bar.m_imag = imag;
        bar.m_real = real;

        writeNew write(bar);
        
        int ID = (int)( 400 - 5*i );
        std::string fileID{ std::to_string(ID) };

        std::string name{ prefix + fileID + suffix };

        H5::H5File file(name, H5F_ACC_TRUNC );
        writeNew writefile(bar);
        writefile.writeMainResults(file);
    }

    selfEnergyValues sev;
    sev.m_at = alphat;
    sev.m_apos = apos;
    sev.m_aneg = aneg;

    H5::H5File file(selfEnergy, H5F_ACC_TRUNC );
        writeNew writefile(sev);
        writefile.writeSe(file);

    std::cout << "File created: \t" << selfEnergy << '\n';
    */

    double tseconds{ t.elapsed() };

    std::cout << "Time taken: " << tseconds << " seconds\n";
    return 0;
}
