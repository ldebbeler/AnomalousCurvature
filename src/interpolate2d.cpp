// source file for 2d Interpolation

#include "interpolate2d.h"

// make sure fxy is correctly ordered in terms of x,y
// seems to be switched, at least for bubble
interpolater2d::interpolater2d(std::vector<double> xa, std::vector<double> ya, std::vector<double> fxy): m_xa(xa), m_ya(ya), m_fxy(fxy)
{
    const gsl_interp2d_type *T = gsl_interp2d_bicubic;  // specify 2dInterpolation type: here bicubic
    const size_t nx = m_xa.size();                        // x grid points 
    const size_t ny = m_ya.size();                        // y grid points
    double *za = new double[nx * ny];                   // allocate memory for known function values on grid

    spline = gsl_spline2d_alloc(T, nx, ny);// create spline object for interpolation, pointer

    xacc = gsl_interp_accel_alloc();    // alternatively, define these in evaluate function
    yacc = gsl_interp_accel_alloc();

    /* set z grid values */
    for(std::size_t i=0; i<nx; i++){
        for(std::size_t j=0; j<ny; j++){
            gsl_spline2d_set(spline, za, i, j, m_fxy[nx*j + i]);
        }
    }

    /* initialize interpolation, this function requires C-style arrays for the grid. Turn vector into this by &vector[0] */
    gsl_spline2d_init(spline, &m_xa[0], &m_ya[0], za, nx, ny);
    delete(za);
}

interpolater2d::~interpolater2d() noexcept{
    gsl_spline2d_free(spline);
    gsl_interp_accel_free(xacc);
    gsl_interp_accel_free(yacc);
}

double interpolater2d::evaluate(double x, double y){
// if x or y out of range, manually interpolate linear or constant
    //gsl_interp_accel *xacc = gsl_interp_accel_alloc();   // state variables for interpolation lookups ???
    //gsl_interp_accel *yacc = gsl_interp_accel_alloc();

    double result;
    int error = gsl_spline2d_eval_e(spline, x, y, xacc, yacc, &result);
    // throw exception if error nonzero
    if(error){
        std::cout << "Error!\n";
    }

    //gsl_interp_accel_free(xacc);
    //gsl_interp_accel_free(yacc);

    return result; 
}


