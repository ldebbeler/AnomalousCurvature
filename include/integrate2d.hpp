# pragma once

# include <vector>
# include <complex>
# include "cubature.h"
# include "d2.hpp"

class integrate2Dmeasure2pi {

	typedef std::complex<double> complex_type;

	// Limit of the BZ
	double kx_min { -M_PI };
	double kx_max {  M_PI };
	double ky_min { -M_PI };
	double ky_max {  M_PI };

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
		*fval = (*func)( double2(k[0],k[1]) ).real();
		return 0;
	}

	template <typename func_t>
	static int integrandIm(unsigned ndim, const double *k, void *param, unsigned fdim, double *fval)
	{
		func_t *func = ((param_t<func_t>*)param)->func;
		*fval = (*func)( double2(k[0],k[1]) ).imag();
		return 0;
	}


	/* Vector integrand */
	template<typename funcv_t>
	struct paramv_t {
		funcv_t *funcv;
		paramv_t(funcv_t *f) : funcv(f){ }
	};

	//template <typename funcv_t>
	//typename std::enable_if< getArgumentCount( funcv_t ) == 2, static int >::type
	//integrandvRe(unsigned ndim, const double *k, void *param, unsigned fdim, double *fval)
	//{
	//	funcv_t *func = ((paramv_t<funcv_t>*)param)->funcv;
	//	for ( size_t i=0; i < fdim ; ++i )
	//		fval[i] = (*func)( double2(k[0],k[1]), i ).real();
	//	return 0;
	//}

	//template <typename funcv_t>
	//typename std::enable_if< getArgumentCount( funcv_t ) == 2, static int >::type
	//integrandvIm(unsigned ndim, const double *k, void *param, unsigned fdim, double *fval)
	//{
	//	funcv_t *func = ((paramv_t<funcv_t>*)param)->funcv;
	//	for ( size_t i=0; i < fdim ; ++i )
	//		fval[i] = (*func)( double2(k[0],k[1]), i ).imag();
	//	return 0;
	//}

	template <typename funcv_t>
	static int integrandvRe(unsigned ndim, const double *k, void *param, unsigned fdim, double *fval)
	{
		funcv_t *func = ((paramv_t<funcv_t>*)param)->funcv;
		auto vec = (*func)( double2(k[0],k[1] ) );
		for ( size_t i=0; i < fdim ; ++i )
			fval[i] = vec[i].real();
		return 0;
	}

	template <typename funcv_t>
	static int integrandvIm(unsigned ndim, const double *k, void *param, unsigned fdim, double *fval)
	{
		funcv_t *func = ((paramv_t<funcv_t>*)param)->funcv;
		auto vec = (*func)( double2(k[0],k[1] ) );
		for ( size_t i=0; i < fdim ; ++i )
			fval[i] = vec[i].imag();
		return 0;
	}


public:

	integrate2Dmeasure2pi( ) { }

	// Setup the BZ limits
	integrate2Dmeasure2pi(const double kkx_min, const double kkx_max, const double kky_min, const double kky_max)
		: kx_min(kkx_min), kx_max(kkx_max), ky_min(kky_min), ky_max(kky_max) { }

	template <typename func_t>
	std::complex<double> integrate( func_t func , double prec = 1.e-4, int steps = 50000)
	{
		double xmin[2], xmax[2];
		//xmin[0] = xmin[1] = -M_PI;
		//xmax[0] = xmax[1] =  M_PI;
		xmin[0] = kx_min;
		xmax[0] = kx_max;
		xmin[1] = ky_min;
		xmax[1] = ky_max;

		param_t<func_t> param(&func);

		// pointer to function
		integrand funcRe = &integrandRe < func_t >;
		integrand funcIm = &integrandIm < func_t >;

		double resRe, resIm, err;
        // absolute and relative error the same: prec
		hcubature(1, funcRe, (void*)&param, 2, xmin, xmax, steps, prec, prec, ERROR_INDIVIDUAL, &resRe, &err);
		hcubature(1, funcIm, (void*)&param, 2, xmin, xmax, steps, prec, prec, ERROR_INDIVIDUAL, &resIm, &err);
		return std::complex<double>(resRe,resIm)/(4.*M_PI*M_PI);
	}

	template <typename funcv_t>
	std::vector<complex_type> integratev( size_t size, funcv_t funcv, double prec = 1.e-4)
	{
		double xmin[2], xmax[2];
		//xmin[0] = xmin[1] = -M_PI;
		//xmax[0] = xmax[1] =  M_PI;
		xmin[0] = kx_min;
		xmax[0] = kx_max;
		xmin[1] = ky_min;
		xmax[1] = ky_max;

		paramv_t<funcv_t> param(&funcv);

		// pointer to function
		integrand funcvRe = &integrandvRe < funcv_t >;
		integrand funcvIm = &integrandvIm < funcv_t >;

		/* real arrays */
		//size_t size = res.size();
		std::vector< double > resRe( size ) , resIm ( size ), err ( size );
		hcubature( size, funcvRe, (void*)&param, 2, xmin, xmax, 50000, prec, prec, ERROR_INDIVIDUAL, resRe.data(), err.data());
		hcubature( size, funcvIm, (void*)&param, 2, xmin, xmax, 50000, prec, prec, ERROR_INDIVIDUAL, resIm.data(), err.data());

		/* too res */
		std::vector<complex_type> res(size);
		for ( size_t i = 0 ; i < size ; ++i )
			res[i] = complex_type ( resRe[i], resIm[i] ) / (4.*M_PI*M_PI);
		return res;
	}



};
