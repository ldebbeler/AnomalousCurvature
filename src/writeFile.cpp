#include "writeFile.h"
#include <boost/multi_array.hpp>
#include <iostream>
#include "constants.h"

writeNew::writeNew(const scalingValues& scaling): m_scaling(scaling) {}

writeNew::writeNew(const selfEnergyValues& sev): m_sev(sev) {}

writeNew::writeNew() {}

void writeNew::write1dvector(H5::H5File file, const std::string& datasetname, const std::vector<double>& vector){
    auto datatype = H5::PredType::NATIVE_DOUBLE;
    datatype.setOrder( H5T_ORDER_LE );
    hsize_t dimv[1];
    dimv[0] = vector.size();
    H5::DataSpace vspace(1, dimv); 

    H5::DataSet dataset = file.createDataSet(datasetname, datatype, vspace);
    dataset.write(&vector[0], datatype);
}

// this probably does not work for an empty vector
void writeNew::write2dvector(H5::H5File file, const std::string& datasetname, const std::vector<std::vector<double>>& vector){
    // replace 2d vector by boost multi_array
    boost::multi_array<double, 2> array( boost::extents[vector.size()][vector[0].size()] );
    for(std::size_t i=0; i<vector.size(); i++){
        for(std::size_t j=0; j<vector[0].size(); j++){
            array[i][j] = vector[i][j];
        }
    }
    auto datatype = H5::PredType::NATIVE_DOUBLE;
    datatype.setOrder( H5T_ORDER_LE );

    std::vector<hsize_t> dimensions(array.shape(), array.shape() +2);
    H5::DataSpace dataspace(2, dimensions.data()); 

    H5::DataSet dataset = file.createDataSet(datasetname, datatype, dataspace);
    dataset.write(array.data(), datatype);
}

void writeNew::createGroup(H5::H5File file, const std::string& groupname){
    file.createGroup(groupname.c_str());
}

void writeNew::writeMainResults(H5::H5File file){
    H5::FloatType double_type(H5::PredType::NATIVE_DOUBLE);  
    H5::FloatType int_type(H5::PredType::NATIVE_INT);  
    H5::DataSpace att_space(H5S_SCALAR);

// group and data for scaling functions
    H5::Group scaling = file.createGroup("/BubbleScale");

    write1dvector(file, "/BubbleScale/Variable", m_scaling.m_x);
    write1dvector(file, "/BubbleScale/Real", m_scaling.m_real);
    write1dvector(file, "/BubbleScale/Imag", m_scaling.m_imag);
    write1dvector(file, "/BubbleScale/ImagDeriv", m_scaling.m_imagDeriv);
    
    H5::Attribute tangential_exponent = scaling.createAttribute( "alphaT", double_type, att_space );
    tangential_exponent.write( double_type, &m_scaling.m_T);

    H5::Group seRad = file.createGroup("/SelfEnergyRadial");

    write1dvector(file, "/SelfEnergyRadial/arPos", m_scaling.m_arPos);
    write1dvector(file, "/SelfEnergyRadial/arNeg", m_scaling.m_arNeg);
    write1dvector(file, "/SelfEnergyRadial/radInterPos", m_scaling.m_radInterPos);
    write1dvector(file, "/SelfEnergyRadial/radInterNeg", m_scaling.m_radInterNeg);

    write1dvector(file, "/SelfEnergyRadial/kr", m_scaling.m_kr);
    write1dvector(file, "/SelfEnergyRadial/seRadPos", m_scaling.m_seRadPos);
    write1dvector(file, "/SelfEnergyRadial/seRadNeg", m_scaling.m_seRadNeg);
 
    H5::Group se = file.createGroup("/SelfEnergy");

    H5::Attribute apos = se.createAttribute("aPos", double_type, att_space);
    apos.write( double_type, &m_scaling.m_apos);
    H5::Attribute aneg = se.createAttribute("aNeg", double_type, att_space);
    aneg.write( double_type, &m_scaling.m_aneg);
  
    H5::Attribute intM = se.createAttribute("M", double_type, att_space);
    intM.write( double_type, &M);
    H5::Attribute intN = se.createAttribute("N", double_type, att_space);
    intN.write( double_type, &N);
 
    H5::Attribute vFermi = se.createAttribute("v", double_type, att_space);
    vFermi.write( double_type, &vF);
    H5::Attribute bParam = se.createAttribute("b", double_type, att_space);
    bParam.write( double_type, &b);

    H5::Group seTang = file.createGroup("/SelfEnergyTangential");

    write1dvector(file, "/SelfEnergyTangential/ktilde", m_scaling.m_ktilde);
    write1dvector(file, "/SelfEnergyTangential/atPos", m_scaling.m_atPos);
    write1dvector(file, "/SelfEnergyTangential/atNeg", m_scaling.m_atNeg);
    write1dvector(file, "/SelfEnergyTangential/interk", m_scaling.m_interk);
    write1dvector(file, "/SelfEnergyTangential/subtractPos", m_scaling.m_subtractPos);
    write1dvector(file, "/SelfEnergyTangential/subtractNeg", m_scaling.m_subtractNeg);
    write1dvector(file, "/SelfEnergyTangential/freqs", m_scaling.m_tangFreqs);
    write1dvector(file, "/SelfEnergyTangential/freqPos", m_scaling.m_tangFreqPos);
    write1dvector(file, "/SelfEnergyTangential/freqNeg", m_scaling.m_tangFreqNeg);
    write1dvector(file, "/SelfEnergyTangential/kt", m_scaling.m_kt);
    write1dvector(file, "/SelfEnergyTangential/seTang", m_scaling.m_seTang);

    H5::Attribute dpos = seTang.createAttribute("DPos", double_type, att_space);
    dpos.write( double_type, &m_scaling.m_DPos);
 
    H5::Attribute dneg = seTang.createAttribute("DNeg", double_type, att_space);
    dneg.write( double_type, &m_scaling.m_DNeg);
 
    H5::Attribute cpos = seTang.createAttribute("CPos", double_type, att_space);
    cpos.write( double_type, &m_scaling.m_CPos);
 
    H5::Attribute cneg = seTang.createAttribute("CNeg", double_type, att_space);
    cneg.write( double_type, &m_scaling.m_CNeg);

    H5::Attribute deltab = seTang.createAttribute("deltab", double_type, att_space);
    deltab.write( double_type, &m_scaling.m_deltab);
 
}

