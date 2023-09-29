// source file for 1d Interpolation
#include <interpolate1d.h>
#include <iostream>


interpolater1d::interpolater1d(std::vector<double> xa, std::vector<double> fx): m_xa(xa), m_fx(fx)
{
    const size_t nx = m_xa.size();                        // x grid points 

    spline = gsl_spline_alloc(gsl_interp_cspline, nx);  // create spline object for interpolation

    // initialize interpolation, this function requires C-style arrays for the grid. Turn vector into this by &vector[0]
    gsl_spline_init(spline, &m_xa[0], &m_fx[0], nx);

    xacc = gsl_interp_accel_alloc();
}

interpolater1d::~interpolater1d() noexcept{
    if(spline==NULL){
        gsl_spline_free(spline);
        gsl_interp_accel_free(xacc);
    }
}

interpolater1d::interpolater1d(interpolater1d& rhs){
    m_xa = rhs.m_xa;
    m_fx = rhs.m_fx;
    //spline = rhs.spline;
    //xacc = rhs.xacc; 
    //rhs.spline = NULL;
    //rhs.xacc = NULL;
    spline = gsl_spline_alloc(gsl_interp_cspline, m_xa.size());  // create spline object for interpolation
    gsl_spline_init(spline, &m_xa[0], &m_fx[0], m_xa.size());
    xacc = gsl_interp_accel_alloc();
 
}

// pass by const ref not possible. Need to change rhs for desired behavior
//interpolater1d::interpolater1d(const interpolater1d& rhs);

double interpolater1d::evaluate(double x){
    //gsl_interp_accel *xacc = gsl_interp_accel_alloc();   // state variables for interpolation lookups 
    double result;

    int error = gsl_spline_eval_e(spline, x, xacc, &result);

    // throw exception if error nonzero
    if (error){         // any nonzero int evaluates to true
        std::cout << "At argument = \t" << x << '\n';
        std::cout << "Could not evaluate interpolated spline." << '\n';
        std::abort();
        //throw interpolatefail();
    }

    return result; 
}

double interpolater1d::firstDevaluate(double x){
    return gsl_spline_eval_deriv(spline, x, xacc);
}

double interpolater1d::secondDevaluate(double x){
    return gsl_spline_eval_deriv2(spline, x, xacc);
}
