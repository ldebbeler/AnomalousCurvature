#ifndef DATASTRUCTS_H
#define DATASTRUCTS_H
#include <vector>

struct scalingValues{
    double m_T;

    // bubble scaling function
    std::vector<double> m_x;       
    std::vector<double> m_real;   
    std::vector<double> m_imag;  
    std::vector<double> m_imagDeriv;

    // radial self energy scaling function
    std::vector<double> m_arPos;
    std::vector<double> m_arNeg;
    std::vector<double> m_radInterPos;
    std::vector<double> m_radInterNeg;

    // radial self energy
    std::vector<double> m_kr;
    std::vector<double> m_seRadPos;
    std::vector<double> m_seRadNeg;

    // tangential self energy scaling function
    std::vector<double> m_ktilde;
    std::vector<double> m_atPos;
    std::vector<double> m_atNeg;
    std::vector<double> m_interk;
    std::vector<double> m_subtractPos;
    std::vector<double> m_subtractNeg;
    std::vector<double> m_tangFreqs;
    std::vector<double> m_tangFreqPos;
    std::vector<double> m_tangFreqNeg;
    double m_DPos;
    double m_DNeg;
    double m_CPos;
    double m_CNeg;
    double m_deltab;

    // tangential self energy
    std::vector<double> m_kt;
    std::vector<double> m_seTang;

    // frequency dependence self energy
    double m_apos;
    double m_aneg;
}; 

struct selfEnergyValues{
    std::vector<double> m_at;
    std::vector<double> m_apos;
    std::vector<double> m_aneg;
};


#endif //DATASTRUCTS_H
