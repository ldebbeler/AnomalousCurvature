//File: constants.h
#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <cmath>
#include <complex>

//const std::string fileName{ "data/foo108.h5" };
const std::string prefix{ "data/expandedTang" };
const std::string suffix{ ".h5" };

//const std::string name{ "data/numBubble3.h5" };
//const std::string selfEnergy{ "data/SelfEnergy.h5" };

//const std::string fileName{ "data/genericOneLoop/subtractZero35.h5" };

inline constexpr std::complex<double> I{ std::complex<double>(0.0,1.0) };   // imaginary unit
inline constexpr double vF{ 1.0 };  // Fermi Velocity
inline constexpr double b{ 1.0 };   // tangential dispersion coefficient
inline constexpr double M{ 3.0 };   // density wave components, 1 for charge 3 for spin
inline constexpr double N{ 2.0 };   // spin multiplicity

// scaling function requires rootfinding for imaginary part and 1d integral for real part
inline constexpr double x_down{ 1e-6 };
inline constexpr double x_up{ 100 };
inline constexpr double epsrel { 1e-3 };
inline constexpr int max_iter{ 100 };
inline constexpr double eps{ 1e-4 };
// 1d integral:
inline constexpr double IR{ 1e-5 }; // lower cutoff for k Integral Bubble (scaling Function and extrapolation coefficients)
inline constexpr double UV{ 1e5 };
inline constexpr int steps{ 200000 };
inline constexpr double precision{ 1e-5 };
// kramers-kronig integral for tangential self energy. Modified cutoff due to subtraction of quadratic term
inline constexpr double IRTang{ 1e-10 };    //fix 1e-10 for now
inline constexpr double UVTang{ 1e8 };      //fix 1e8

// assembly of scaling function. Interpolation in [-Lambda, Lambda], extrapolation beyond
inline constexpr int nodes{ 5001 };
//inline constexpr int nodes{ 51 };
inline constexpr double Lambda{ 500.0 };        // 100 leads to smooth tangential scaling function
// bubble calculation via scaling function
// specification for 2d integral (hcubature)
inline constexpr double prec2d{ 1e-10 };
inline constexpr int steps2d{ 5000*1000 };
inline constexpr double IR2d{ 1e-7 };
inline constexpr double UV2d{ 1e5 };

// exact quartic scaling function used for bubble
inline constexpr bool exactQuartic{ false };
// scaling function for radial self energy
inline constexpr int nodesR{ 201 };
inline constexpr double LambdaR{ 5.0 };

// calibrate with these values + range for quartic fit in newSe
inline constexpr int nodesT{ 301 };         // 201 seems good
inline constexpr double LambdaT{ 5.0 };     // 5.0 seems good, 7.0 as well
// regression for quartic fit
inline constexpr int quarticNodes{ 200 };   // 200 seems like a good value
inline constexpr double maxValue{ 1.0 };    // careful, this is very sensitive. 1.0
// value below which quartic approx
inline constexpr double quart{ 0.5 };       // 0.7, below alpha=2.5 should be decreased e.g. 0.5

/*
// specification for grid of scaling funcion
inline constexpr bool logSpacing{ true }; //grid points of scaling functions have logarithmic spacing. Otherwise linear spacing
inline constexpr std::size_t N{ 30 };    //Number of grid points for scaling functions
inline constexpr double minx{ 1e-2 };     //smallest nonzero grid value for logarithmic spacing
inline constexpr double extraPolate{ 100.0 }; // scaling function is extrapolated for values larger than this
*/

// self energy calculation imaginary part via 2d q-integration (hcubature)
inline constexpr int steps2dq{ 1000*500 };
inline constexpr double prec2dq{ 1e-12 };
inline constexpr double IRq{ 1e-8 };
inline constexpr double UVq{ 1e5 };

#endif  //CONSTANTS_H
