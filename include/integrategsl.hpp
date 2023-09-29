#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

#include "constants.h"

class integrater{

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
    integrater() {}

    template<typename func_t>
    // a: lower boundary, b: upper boundary
    double integrate( func_t func, double a, double b)
    {
        gsl_integration_workspace * w = gsl_integration_workspace_alloc (steps);

        double result;
        double error;

        gsl_function_pp<decltype(func)> Fp(func);           // construct an object of class gsl_function_pp (that specific class with the right func type of the template)
        gsl_function *F = static_cast<gsl_function*>(&Fp);

        int status = gsl_integration_qags (F, a, b, precision, precision, steps, w, &result, &error);
        gsl_integration_workspace_free (w);

        // throw exception here if return value nonzero
        if(status){ throw 42;  }
        return result;
    }
};

