//File: writeNew.h
#ifndef WRITENEW_H
#define WRITENEW_H
#include <H5Cpp.h>
#include <vector>
#include <string>
#include "datastructs.h"

class writeNew {
    scalingValues m_scaling;
    coefficients m_sev;

public:
    writeNew(const scalingValues& scaling);
    writeNew(const coefficients& sev);
    writeNew();

private:
    void write1dvector(H5::H5File file, const std::string& datasetname, const std::vector<double>& vector);  // create dataset with 1d vector

    void write2dvector(H5::H5File file, const std::string& datasetname, const std::vector<std::vector<double>>& vector); // create dataset with 2d vector

    void createGroup(H5::H5File file, const std::string& groupname);    

public:
    void writeMainResults(H5::H5File file);  

    void writeSe(H5::H5File file);

    void writeCoefficients(H5::H5File file);
};    

#endif //WRITENEW_H
