#ifndef WRITEH5_H
#define WRITEH5_H
#include <H5Cpp.h>
#include <vector>
#include <string>
#include "datastructs.h"

class writeh5 {
    scalingValues m_scaling;
    bubbleValues m_bubble;
    selfEnergyValues m_selfEnergy;

public:
    writeh5(const scalingValues& scaling, const bubbleValues& bubble, const selfEnergyValues& sev);
    writeh5(const scalingValues& scaling, const selfEnergyValues& sev);
    writeh5();

private:
    void write1dvector(H5::H5File file, const std::string& datasetname, const std::vector<double>& vector);  // create dataset with 1d vector

    void write2dvector(H5::H5File file, const std::string& datasetname, const std::vector<std::vector<double>>& vector); // create dataset with 2d vector

    void createGroup(H5::H5File file, const std::string& groupname);    // create group with dataset m_cutoffs, and subgroups "main results"and "flow"

public:
    void writeMainResults(H5::H5File file);    //write most important results into this group: selfenergy and static susceptibility at end of flow. Z-factors during flow.

    void writeFunc(H5::H5File file, const std::vector<double>& kr, const std::vector<double>& funcPos, const std::vector<double>& funcNeg);
};    

#endif //WRITEH5_H
