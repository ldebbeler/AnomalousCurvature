#ifndef INTERPOLATE2D_H
#define INTERPOLATE2D_H
// Header file for GSL-based 2d Interpolation

#include <vector>
#include <iostream>

#include <gsl/gsl_math.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

//#include "exceptclass.h"

class interpolater2d{
    std::vector<double> m_xa; // 1st grid of values
    std::vector<double> m_ya; // 2nd grid of values
    std::vector<double> m_fxy;// function values on underlying grid
    gsl_spline2d* spline;     // pointer to spline object. Necessary as interpolation function require pointer as argument.
    gsl_interp_accel* xacc;   // alternatively create and destroy these locally in evaluate function
    gsl_interp_accel* yacc;

public: 
    // construct spline, that interpolates between given values
    interpolater2d(std::vector<double> xa, std::vector<double> ya, std::vector<double> fxy);

    // destructor. Clean up spline memory
    ~interpolater2d() ;

    // evaluate spline at any 2d point
    double evaluate(double x, double y);
};

#endif
