#ifndef INTERPOLATE1D_H
#define INTERPOLATE1D_H
// Header File for GSL-based 1D Interpolation

#include <vector>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>

class interpolater1d{
    std::vector<double> m_xa; // 1st grid of values
    std::vector<double> m_fx; // function values on underlying grid
    gsl_spline* spline;       // pointer to spline object. Necessary as interpolation function require pointer as argument.
    gsl_interp_accel* xacc;

public: 
    interpolater1d(std::vector<double> xa, std::vector<double> fx);

    ~interpolater1d() ;

    interpolater1d(interpolater1d& rhs);
    // we cannot pass by const ref for copy assignment
    //interpolater1d(const interpolater1d& rhs);
private:
    // we do not want to copy objects of this class. Therefore we declare these methods private and dont define them
    // this is actually not true. Copying is complicated because pointers unique?
    interpolater1d& operator=(const interpolater1d&);

public:
    double evaluate(double x);
    double firstDevaluate(double x);
    double secondDevaluate(double x);
};

#endif
