# pragma once

# include <vector>
# include <complex>
# include <tuple>
# include <iostream>
# include "cubature.h"
# include "dcuhre_interface.hpp"
# include "d2.hpp"

typedef enum {
    lib_hcubature,
    lib_dcuhre
} integrallib_type;

class averageBZ {

	typedef std::complex<double> complex_type;

	// Limit of the BZ
	double kx_min { -M_PI };
	double kx_max {  M_PI };
	double ky_min { -M_PI };
	double ky_max {  M_PI };

	template<typename func_t>
	struct param_t {
                double theta { 0. };
		func_t *func;
		param_t(func_t *f) : func(f){ }
	};

	template <typename R, typename ... Types> constexpr size_t getArgumentCount( R(*f)(Types ...))
	{ return sizeof...(Types); }

	template <typename func_t>
	static void integrand_dcuhre_Re(double *k, int ndim, double *fval, int fdim, void *param)
	{
	    func_t *func = ((param_t<func_t>*)param)->func;
	    *fval = (*func)( double2(k[0], k[1]) ).real();
	}

	template <typename func_t>
	static void integrand_dcuhre_Im(double *k, int ndim, double *fval, int fdim, void *param)
	{
	    func_t *func = ((param_t<func_t>*)param)->func;
	    *fval = (*func)( double2(k[0], k[1]) ).imag();
	}


	template <typename func_t>
	static int integrand_hcubature_Re(unsigned ndim, const double *k, void *param, unsigned fdim, double *fval)
	{
		func_t *func = ((param_t<func_t>*)param)->func;
		*fval = (*func)( double2(k[0],k[1]) ).real();
		return 0;
	}

	template <typename func_t>
	static int integrand_hcubature_Im(unsigned ndim, const double *k, void *param, unsigned fdim, double *fval)
	{
		func_t *func = ((param_t<func_t>*)param)->func;
		*fval = (*func)( double2(k[0],k[1]) ).imag();
		return 0;
	}

        template <typename func_t>
        static int integrand_radial_Re(unsigned ndim, const double *r, void *param, unsigned fdim, double *fval)
        {
            *fval = integrand_radial<func_t>(ndim, r, param, fdim).real();
            return 0;
        }

        template <typename func_t>
        static int integrand_radial_Im(unsigned ndim, const double *r, void *param, unsigned fdim, double *fval)
        { 
            *fval = integrand_radial<func_t>(ndim, r, param, fdim).imag();
            return 0;
        }

        template <typename func_t>
        static complex_type integrand_radial(unsigned ndim, const double *r, void *param, unsigned fdim)
        {
                param_t<func_t>* p = (param_t<func_t>*)param;

                // rad 
                double theta = p->theta;
                func_t *func = p->func;

                // get cartesian coordinate
                double kx = r[0] * cos(theta);
                double ky = r[0] * sin(theta);
                double Jac= r[0];

                return Jac * (*func)( double2(kx, ky) );
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

        /* INTEGRATE A MATRIX */
	/* Vector integrand */
	template<typename funcv_t>
	struct paramm_t {
		funcv_t *funcv;
                size_t mat_size;
		paramm_t(funcv_t *f, size_t ms) : funcv(f), mat_size(ms){ }
	};



	/* integrate a double */
	template <typename func_t>
	static int integrand_double(unsigned ndim, const double *k, void *param, unsigned fdim, double *fval)
	{
		func_t *func = ((param_t<func_t>*)param)->func;
		*fval = (*func)( double2(k[0],k[1]) );
		return 0;
	}

        // In case of radial integration: we need (for fixed theta) the extrema for r variable
        auto get_radial_extrema(double theta_in)
        {
            // 4 quadranti
            //
            //  |---------------|
            //  |       2       |
            //  |               |
            //  |               |
            //  | 3           1 |
            //  |               |
            //  |               |
            //  |       4       |
            //  |---------------|
            auto [ ax1, ax2, ay1, ay2 ] = std::make_tuple( std::abs(kx_min), std::abs(kx_max), std::abs(ky_min), std::abs(ky_max) );

            double theta = fmod( theta_in, 2.*M_PI );

            // Get limit of first quadrant
            double theta_1 =  atan( ay2 / ax2 );
            double theta_4 = 2.*M_PI - atan( ay1 / ax2 );
            double theta_2 = M_PI - atan( ay2 / ax1 );
            double theta_3 = M_PI + atan( ay1 / ax1 );

            if( theta <= theta_1 or (theta > theta_4 and theta < 2.*M_PI + theta_1)  ) {

                double loc_theta = theta;
                // I am in the first quadrant
                double rmax = sqrt( ax2*ax2 * ( 1. + tan(loc_theta)*tan(loc_theta) ) );
                //std::cout << "Rma1 " << theta << "\t" << rmax << std::endl;
                return std::make_tuple( 0., rmax );

            }
            else if( theta >  theta_1 and theta <= theta_2 ) {

                double loc_theta = M_PI/2. - theta;

                // I am in the second quadrant
                double rmax = sqrt( ay2*ay2 * ( 1. + tan(loc_theta)*tan(loc_theta) ) );
                //std::cout << "Rma2 " << theta << "\t" << rmax << std::endl;
                return std::make_tuple( 0., rmax );

            }
            else if( theta >  theta_2 and theta <= theta_3 ) {

                double loc_theta = theta - M_PI;

                // I am in the third quadrant
                double rmax = sqrt( ax1*ax1 * ( 1. + tan(loc_theta)*tan(loc_theta) ) );
                //std::cout << "Rma3 " << theta << "\t" << rmax << std::endl;
                return std::make_tuple( 0., rmax );

            }
            else if( theta >  theta_3 and theta <= theta_4 ) {

                double loc_theta = theta - 3./2.*M_PI;

                // I am in the third quadrant
                double rmax = sqrt( ay1*ay1 * ( 1. + tan(loc_theta)*tan(loc_theta) ) );
                //std::cout << "Rma4 " << theta << "\t" << loc_theta << "\t" << rmax << std::endl;
                return std::make_tuple( 0., rmax );
            }
            else {
                std::cerr << ( "Why here? get_radial_extrema internal error" ) << std::endl;
                return std::make_tuple( 0., 0. );
            }
            
        }


public:

	averageBZ( ) { }

	// Setup the BZ limits
	averageBZ(const double kkx_min, const double kkx_max, const double kky_min, const double kky_max)
		: kx_min(kkx_min), kx_max(kkx_max), ky_min(kky_min), ky_max(kky_max) { }

        template <typename func_t>
        std::complex<double> integrate( func_t func, double prec = 1.e-4 , int maxEv=5e4)
        {
            return integrate_hcubature(func, prec, maxEv);
        }

        template <typename func_t>
        std::complex<double> integrate( func_t func, integrallib_type lib, double prec = 1.e-4, int maxEv = 5e4 )
        {
            if ( lib == lib_hcubature )
                return integrate_hcubature(func, prec, maxEv);
            else if( lib == lib_dcuhre )
                return integrate_dcuhre(func, prec);
            else {
                std::cerr << ("Unrecognized integration library. Possibilities: hcubature, dcuhre") << std::endl;
                std::abort();                
                return std::complex<double>( std::nan("1"), std::nan("1") );
            }
        }

	template <typename func_t>
	std::complex<double> integrate_dcuhre( func_t func , double prec = 1.e-4)
	{
		//auto vol = (kx_max - kx_min) * (ky_max - ky_min);
        auto vol = 2*M_PI*2*M_PI;

		param_t<func_t> param(&func);

		// pointer to function
		integrand_dcuhre funcRe = &integrand_dcuhre_Re < func_t >;
		integrand_dcuhre funcIm = &integrand_dcuhre_Im < func_t >;

		double resRe[1], resIm[1];
                int errRe, errIm;

		dcuhre_interface_2D(1, kx_min, kx_max, ky_min, ky_max, prec, prec, (double*)resRe, funcRe, (void*)&param, &errRe); 
		dcuhre_interface_2D(1, kx_min, kx_max, ky_min, ky_max, prec, prec, (double*)resIm, funcIm, (void*)&param, &errIm); 

                if ( errRe != 0 )
                    std::cerr << "dcuhre Real gave the return value: " << errRe << std::endl;
                if ( errIm != 0 )
                    std::cerr << "dcuhre Real gave the return value: " << errIm << std::endl;

		return std::complex<double>(resRe[0],resIm[0])/vol;
	}



	template <typename func_t>
	std::complex<double> integrate_hcubature( func_t func , double prec = 1.e-4, int maxEv=5e4)
	{
		double xmin[2], xmax[2];
		xmin[0] = kx_min;
		xmax[0] = kx_max;
		xmin[1] = ky_min;
		xmax[1] = ky_max;
		auto vol = (kx_max - kx_min) * (ky_max - ky_min);

		param_t<func_t> param(&func);

		// pointer to function
		integrand funcRe = &integrand_hcubature_Re < func_t >;
		integrand funcIm = &integrand_hcubature_Im < func_t >;

		double resRe, resIm, errRe, errIm;

		hcubature(1, funcRe, (void*)&param, 2, xmin, xmax, maxEv, prec, prec, ERROR_INDIVIDUAL, &resRe, &errRe);
		hcubature(1, funcIm, (void*)&param, 2, xmin, xmax, maxEv, prec, prec, ERROR_INDIVIDUAL, &resIm, &errIm);

                //if (errRe > prec) std::cerr << "integrate: error real larger than precision!" << errRe << "\t" << prec << std::endl;
                //if (errIm > prec) std::cerr << "integrate: error imag larger than precision!" << errIm << "\t" << prec << std::endl;

		return std::complex<double>(resRe,resIm)/vol;
	}

	template <typename funcv_t>
	std::vector<complex_type> integratev( size_t size, funcv_t funcv, double prec = 1.e-4)
	{
		double xmin[2], xmax[2];
		xmin[0] = kx_min;
		xmax[0] = kx_max;
		xmin[1] = ky_min;
		xmax[1] = ky_max;
		auto vol = (kx_max - kx_min) * (ky_max - ky_min);

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
			res[i] = complex_type ( resRe[i], resIm[i] ) / vol;
		return res;
	}

	/* integrate a double */
	template <typename func_t>
	double integrate_double( func_t func , double prec = 1.e-4, int maxEv=5e4)
	{
		double xmin[2], xmax[2];
		xmin[0] = kx_min;
		xmax[0] = kx_max;
		xmin[1] = ky_min;
		xmax[1] = ky_max;
		auto vol = (kx_max - kx_min) * (ky_max - ky_min);

		param_t<func_t> param(&func);

		// pointer to function
		integrand funcInte = &integrand_double < func_t >;

		double resRe, err;
		hcubature(1, funcInte, (void*)&param, 2, xmin, xmax, maxEv, prec, prec, ERROR_INDIVIDUAL, &resRe, &err);
		return resRe/vol;
	}


        template <typename func_t>
        auto integrate_radial( func_t func, int num_theta = 100, double prec = 1.e-4, int maxEv=5e4 )
        {
            //constexpr size_t num_theta = 100;

	    auto vol = (kx_max - kx_min) * (ky_max - ky_min);

            // Define the function pointer
            param_t<func_t> param(&func);
            integrand funcRe = &integrand_radial_Re<func_t>;
            integrand funcIm = &integrand_radial_Im<func_t>;

            // Theta loop
            complex_type res = 0.;
            for(int itheta=0; itheta<num_theta; ++itheta) {

                double theta        = itheta * 2.*M_PI / num_theta;
                double delta_theta  = 2.*M_PI / num_theta;
                param.theta         = theta;

                // Calculate the extrema for r
                auto [ rmin, rmax ] = get_radial_extrema(theta);

                double radial_integralRe = 0., radial_integralIm = 0., errRe = 0., errIm = 0.;
                hcubature(1, funcRe, (void*)&param, 1, &rmin, &rmax, maxEv, prec, prec, ERROR_INDIVIDUAL, &radial_integralRe, &errRe);
                hcubature(1, funcIm, (void*)&param, 1, &rmin, &rmax, maxEv, prec, prec, ERROR_INDIVIDUAL, &radial_integralIm, &errIm);
                //std::cout << "Got rmin, rmax, integral: " << rmin << "\t" << rmax << "\t" << radial_integral << "\t" << 0.5*rmax*rmax << std::endl;

                res += delta_theta * complex_type( radial_integralRe, radial_integralIm );
            }
            return res / vol;
        }

        // Here, we suppose that the integrand is highly symmetric. we perform the integral only onver:
        //
        //
        //
        //         /|
        //        / |
        //       /  |
        //      /   |
        //     /    |
        //    / Qui |
        //   /      |
        //  /       |
        // |---------
        //
        //template <typename func_t>
        //auto integrate_radial_supersymmetric( func_t func, double prec = 1.e-4, int maxEv=5e4 )
        //{
        //    constexpr size_t num_theta   = 1000;
        //    constexpr double Total_angle = M_PI/4.; // Supersymmetric

	//    auto vol = (kx_max - kx_min) * (ky_max - ky_min);

        //    // Define the function pointer
        //    param_t<func_t> param(&func);
        //    integrand funcRe = &integrand_radial_Re<func_t>;
        //    integrand funcIm = &integrand_radial_Im<func_t>;

        //    // Theta loop
        //    complex_type res = 0.;
        //    for(int itheta=0; itheta<num_theta; ++itheta) {

        //        double theta        = itheta * Total_angle / num_theta;
        //        double delta_theta  = Total_angle / num_theta;
        //        param.theta         = theta;

        //        // Calculate the extrema for r
        //        auto [ rmin, rmax ] = get_radial_extrema(theta);

        //        double radial_integralRe = 0., radial_integralIm, errRe = 0., errIm = 0.;
        //        hcubature(1, funcRe, (void*)&param, 1, &rmin, &rmax, maxEv, prec, prec, ERROR_INDIVIDUAL, &radial_integralRe, &errRe);
        //        hcubature(1, funcIm, (void*)&param, 1, &rmin, &rmax, maxEv, prec, prec, ERROR_INDIVIDUAL, &radial_integralIm, &errIm);
        //        //std::cout << "Got rmin, rmax, integral: " << rmin << "\t" << rmax << "\t" << radial_integral << "\t" << 0.5*rmax*rmax << std::endl;

        //        res += delta_theta * complex_type( radial_integralRe, radial_integralIm );
        //    }
        //    return 8.*res / vol;
        //}



};
