# pragma once

# include <vector>
# include <complex>
# include "cubature.h"
# include "d3.hpp"

class integrate3Dmeasure2pi {

	typedef std::complex<double> complex_type;

	// Limit of the BZ
	double kx_min { -M_PI };
	double kx_max {  M_PI };
	double ky_min { -M_PI };
	double ky_max {  M_PI };
    double kz_min { -M_PI };
    double kz_max {  M_PI };

	template<typename func_t>
	struct param_t {
		func_t *func;
		param_t(func_t *f) : func(f){ }
	};


	template <typename R, typename ... Types> constexpr size_t getArgumentCount( R(*f)(Types ...))
	{ return sizeof...(Types); }

	template <typename func_t>
	static int integrandRe(unsigned ndim, const double *k, void *param, unsigned fdim, double *fval)
	{
		func_t *func = ((param_t<func_t>*)param)->func;
		*fval = (*func)( double3(k[0],k[1],k[2]) ).real();
		return 0;
	}

	template <typename func_t>
	static int integrandIm(unsigned ndim, const double *k, void *param, unsigned fdim, double *fval)
	{
		func_t *func = ((param_t<func_t>*)param)->func;
		*fval = (*func)( double3(k[0],k[1],k[2]) ).imag();
		return 0;
	}

public:

	integrate3Dmeasure2pi( ) { }

	// Setup the BZ limits
	integrate3Dmeasure2pi(const double kkx_min, const double kkx_max,
        const double kky_min, const double kky_max, const double kkz_min, const double kkz_max):
		    kx_min(kkx_min), kx_max(kkx_max), ky_min(kky_min), ky_max(kky_max),
            kz_min(kkz_min), kz_max(kkz_max) { }

	template <typename func_t>
	std::complex<double> integrate( func_t func , double prec = 1e-4, int steps=50000)
	{
		double xmin[3], xmax[3];
		xmin[0] = kx_min;
		xmax[0] = kx_max;
		xmin[1] = ky_min;
		xmax[1] = ky_max;
        xmin[2] = kz_min;
        xmax[2] = kz_max;

		param_t<func_t> param(&func);

		// pointer to function
		integrand funcRe = &integrandRe < func_t >;
		integrand funcIm = &integrandIm < func_t >;

		double resRe, resIm, err;
        // absolute and relative error the same: prec
		hcubature(1, funcRe, (void*)&param, 3, xmin, xmax, steps, prec, prec, ERROR_INDIVIDUAL, &resRe, &err);
		hcubature(1, funcIm, (void*)&param, 3, xmin, xmax, steps, prec, prec, ERROR_INDIVIDUAL, &resIm, &err);
		return std::complex<double>(resRe,resIm)/(8.*M_PI*M_PI*M_PI);
	}

};
