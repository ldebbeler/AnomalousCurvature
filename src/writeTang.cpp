#include "writeTang.h"
#include <boost/multi_array.hpp>
#include "constants.h"

writeTang::writeTang(const scalingValues& scaling, const bubbleValues& bubble, const seValuesReal& sev):
    m_scaling(scaling), m_bubble(bubble), m_selfEnergy(sev) {}

void writeTang::write1dvector(H5::H5File file, const std::string& datasetname, const std::vector<double>& vector){
    auto datatype = H5::PredType::NATIVE_DOUBLE;
    datatype.setOrder( H5T_ORDER_LE );
    hsize_t dimv[1];
    dimv[0] = vector.size();
    H5::DataSpace vspace(1, dimv); 

    H5::DataSet dataset = file.createDataSet(datasetname, datatype, vspace);
    dataset.write(&vector[0], datatype);
}

void writeTang::write2dvector(H5::H5File file, const std::string& datasetname, const std::vector<std::vector<double>>& vector){
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

void writeTang::createGroup(H5::H5File file, const std::string& groupname){
    file.createGroup(groupname.c_str());
}

void writeTang::writeMainResults(H5::H5File file){
    H5::FloatType double_type(H5::PredType::NATIVE_DOUBLE);  
    H5::FloatType int_type(H5::PredType::NATIVE_INT);  
    H5::DataSpace att_space(H5S_SCALAR);

// group and data for scaling functions
    H5::Group scaling = file.createGroup("/ScalingFunction");

    write1dvector(file, "/ScalingFunction/Variable", m_scaling.m_x);
    write1dvector(file, "/ScalingFunction/Real", m_scaling.m_real);
    write1dvector(file, "/ScalingFunction/Imag", m_scaling.m_imag);

    H5::Attribute upper_cutoff = scaling.createAttribute( "Upper Cutoff Bubble", double_type, att_space );
    upper_cutoff.write( double_type, &UVCutoff);

    H5::Attribute lower_cutoff = scaling.createAttribute( "Lower Cutoff Bubble", double_type, att_space );
    lower_cutoff.write( double_type, &IRCutoff);

    H5::Attribute precision = scaling.createAttribute( "Precision", double_type, att_space );
    precision.write( double_type, &precision);

    H5::Attribute steps = scaling.createAttribute( "Steps", int_type, att_space );
    steps.write( int_type, &steps);

    H5::Attribute vF_const = scaling.createAttribute( "vF", double_type, att_space );
    vF_const.write( double_type, &vF);

    H5::Attribute b_const = scaling.createAttribute( "b", double_type, att_space );
    b_const.write( double_type, &b);

    H5::Attribute extrap = scaling.createAttribute( "extraPolate", double_type, att_space );
    extrap.write( double_type, &extraPolate);

// group and data for bubble with real frequencies
    H5::Group bubbleTest = file.createGroup("/Bubble");
    write1dvector(file, "/Bubble/omega", m_bubble.m_omega);
    write1dvector(file, "/Bubble/omegReal", m_bubble.m_omegReal);
    write1dvector(file, "/Bubble/omegImag", m_bubble.m_omegImag);

    write1dvector(file, "/Bubble/qr", m_bubble.m_qr);
    write1dvector(file, "/Bubble/radReal", m_bubble.m_radReal);
    write1dvector(file, "/Bubble/radImag", m_bubble.m_radImag);
    
    write1dvector(file, "/Bubble/qt", m_bubble.m_qt);
    write1dvector(file, "/Bubble/tangReal", m_bubble.m_tangReal);
    write1dvector(file, "/Bubble/tangImag", m_bubble.m_tangImag);

    H5::Attribute omegqr = bubbleTest.createAttribute( "Freq qr", double_type, att_space );
    omegqr.write( double_type, &m_bubble.m_omegqr);
    H5::Attribute omegqt = bubbleTest.createAttribute( "Freq qt", double_type, att_space );
    omegqt.write( double_type, &m_bubble.m_omegqt);

    H5::Attribute radOmega = bubbleTest.createAttribute( "Radial omega", double_type, att_space );
    radOmega.write( double_type, &m_bubble.m_radOmega);
    H5::Attribute radqt = bubbleTest.createAttribute( "Radial qt", double_type, att_space );
    radqt.write( double_type, &m_bubble.m_radqt);

    H5::Attribute tangOmega = bubbleTest.createAttribute( "Tangential omega", double_type, att_space );
    tangOmega.write( double_type, &m_bubble.m_tangOmega);
    H5::Attribute tangqr = bubbleTest.createAttribute( "Tangential qr", double_type, att_space );
    tangqr.write( double_type, &m_bubble.m_tangqr);

    H5::Group selfEnergy = file.createGroup("/SelfEnergy");
    write1dvector(file, "/SelfEnergy/cutoffs", m_selfEnergy.m_val);
    write1dvector(file, "/SelfEnergy/Ztang", m_selfEnergy.m_SEPos);

}

