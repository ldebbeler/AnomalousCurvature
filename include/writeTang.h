#ifndef WRITETANG_H
#define WRITETANG_H
#include <H5Cpp.h>
#include <vector>
#include <string>
#include "datastructs.h"

class writeTang {
    scalingValues m_scaling;
    bubbleValues m_bubble;
    seValuesReal m_selfEnergy;

public:
    writeTang(const scalingValues& scaling, const bubbleValues& bubble, const seValuesReal& sev);

private:
    void write1dvector(H5::H5File file, const std::string& datasetname, const std::vector<double>& vector);  // create dataset with 1d vector

    void write2dvector(H5::H5File file, const std::string& datasetname, const std::vector<std::vector<double>>& vector); // create dataset with 2d vector

    void createGroup(H5::H5File file, const std::string& groupname);    // create group with dataset m_cutoffs, and subgroups "main results"and "flow"

public:
    void writeMainResults(H5::H5File file);    //write most important results into this group: selfenergy and static susceptibility at end of flow. Z-factors during flow.
};    

#endif //WRITETANG_H
