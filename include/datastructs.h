#ifndef DATASTRUCTS_H
#define DATASTRUCTS_H
#include <vector>

struct scalingValues{
    double m_T;

    // bubble scaling function
    std::vector<double> m_x;
    std::vector<double> m_real;
    std::vector<double> m_imag;

    // frequency dependence self energy
    double m_apos;
    double m_aneg;

    // radial self energy scaling function
    std::vector<double> m_krt;
    std::vector<double> m_arPos;
    std::vector<double> m_arNeg;

    // real part radial self energy
    std::vector<double> m_kr;
    std::vector<double> m_seRad;

    // tangential self energy scaling function
    std::vector<double> m_ktt;
    std::vector<double> m_atPos;
    std::vector<double> m_atNeg;

    //real part tangential self energy
    std::vector<double> m_kt;
    std::vector<double> m_seTang;

    double m_deltab;
    double m_Cpos;
    double m_Cneg;
    double m_Dpos;
    double m_Dneg;

    /*
    std::vector<double> m_interk;
    std::vector<double> m_subtractPos;
    std::vector<double> m_subtractNeg;
    std::vector<double> m_tangFreqs;
    std::vector<double> m_tangFreqPos;
    std::vector<double> m_tangFreqNeg;
    */

}; 

struct selfEnergyValues{
    std::vector<double> m_at;
    std::vector<double> m_apos;
    std::vector<double> m_aneg;
};


#endif //DATASTRUCTS_H
