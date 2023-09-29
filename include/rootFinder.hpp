#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
#include <iostream>

#include "constants.h"

class rootFinder{

    // nested class to bring any function into required gsl_function form
    template< typename F >
        class gsl_function_pp : public gsl_function {       // nested class inheriting from base class (gsl_function)
            public:
            gsl_function_pp(const F& func) : _func(func) {  // Constructor nested class
                                                            // data type gsl_function contains "double (*function) (double x, void * params)" and "void * params"
                function = &gsl_function_pp::invoke;        // invoke static function, can be called without instance of class. Defined below 
                params = this;                              // adress of this object
            }                                                                   
            private:                                                            
            const F& _func;                                         // nested-class member
            static double invoke(double x, void * params) {         // nested-class member function. Static as it does not refer to any instance members
                return static_cast<gsl_function_pp*>(params)->_func(x);         
            }                                                                   
        };       
public: 
    rootFinder() {}     // constructor

    template<typename func_t>
    // x_lo: lower boundary, x_hi: upper boundary, r: start for root
    double findRoot( func_t func, double x_lo, double x_hi, double r)
    {
        const gsl_root_fsolver_type* T{ gsl_root_fsolver_brent };  // specify root solving algorithm
        gsl_root_fsolver* s{ gsl_root_fsolver_alloc(T) };          // solver state for algorithm
        
        // construct an object of class gsl_function_pp (that specific class with the right func type of the template)
        gsl_function_pp<decltype(func)> Fp(func);                   
        gsl_function *F = static_cast<gsl_function*>(&Fp);

        gsl_root_fsolver_set (s, F, x_lo, x_hi);                    // initialize solver to work with function F

        int status;
        int iter{ 0 };

        do          // loop for iteration
        {
            iter++;
            status = gsl_root_fsolver_iterate(s);       // perform iteration of solver
            r = gsl_root_fsolver_root(s);               // current estimate for root
            x_lo = gsl_root_fsolver_x_lower(s);         // current bracketing interval lower end
            x_hi = gsl_root_fsolver_x_upper(s);         // current bracketing interval upper end
            status = gsl_root_test_interval(x_lo, x_hi, 0, epsrel); // tests convergence. GSL_SUCCESS returned when converged
                                                                   // GSL_CONTINUE otherwise
            if(status == GSL_SUCCESS){
                //std::cout << "Converged \n";
            }
        }
        while (status ==GSL_CONTINUE && iter < max_iter);

        gsl_root_fsolver_free (s);

        return r;
    }

};

/*
class straddleCheck{ // check whether there is a sign change in provided interval for rootFinder

    // nested class to bring any function into required gsl_function form
    template< typename F >
        class gsl_function_pp : public gsl_function {       // nested class inheriting from base class (gsl_function)
            public:
            gsl_function_pp(const F& func) : _func(func) {  // Constructor nested class
                                                            // data type gsl_function contains "double (*function) (double x, void * params)" and "void * params"
                function = &gsl_function_pp::invoke;        // invoke static function, can be called without instance of class. Defined below 
                params = this;                              // adress of this object
            }                                                                   
            private:                                                            
            const F& _func;                                         // nested-class member
            static double invoke(double x, void * params) {         // nested-class member function. Static as it does not refer to any instance members
                return static_cast<gsl_function_pp*>(params)->_func(x);         
            }                                                                   
        };       
public: 
    straddleCheck() {}     // constructor

    template<typename func_t>
    // x_lo: lower boundary, x_hi: upper boundary, r: start for root
    bool signChange( func_t func, double x_lo, double x_hi, double r)
    {
        gsl_function_pp<decltype(func)> Fp(func);           // construct an object of class gsl_function_pp (that specific class with the right func type of the template)
        gsl_function *F = static_cast<gsl_function*>(&Fp);

        int status;

        return true;
    }

};
*/
