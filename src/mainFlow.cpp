#include <iostream>
#include "scaling.h"
#include "grid.h"
#include "bubble.h"
#include "interpolate1d.h"
#include "datastructs.h"
#include "testBubble.h"
#include "writeTang.h"
#include "tangFlow.h"
#include "timer.hpp"
#include "constants.h"
#include <boost/numeric/odeint.hpp>

struct push_back_state_and_time
{
    std::vector<std::vector<double>>& m_states;
    std::vector<double>& m_times;

    push_back_state_and_time(std::vector<std::vector<double>>& states, std::vector<double>& times):
        m_states(states), m_times(times) {}
    
    void operator()( const std::vector<double>& x, double l)
    {
        m_states.push_back( x );
        m_times.push_back(std::exp(l));
    }
};

int main() {

    // save bubble values in these structs
    scalingValues scale;
    bubbleValues bv; 
    bubble bub;
    testBubble foo(bub);
    foo.writeBubble(scale, bv);

    tangFlow bar(bub);

    Timer t;
    std::cout << bar.rightFlow(0.1) << '\n';
    std::cout << bar.rightFlow(0.2) << '\n';
    std::cout << bar.rightFlow(0.3) << '\n';
    std::cout << bar.rightFlow(-0.1) << '\n';
    std::cout << bar.rightFlow(-0.2) << '\n';
    std::cout << bar.rightFlow(-0.3) << '\n';

    using namespace boost::numeric::odeint;

    std::vector<double> x(1);
    x[0] = 0.0;

    std::vector<std::vector<double>> x_vec;
    std::vector<double> cutoffs;

    typedef runge_kutta_cash_karp54< std::vector<double> > error_stepper_type;

    typedef controlled_runge_kutta< error_stepper_type > controlled_stepper_type;

    double abs_err = 1.0e-3 , rel_err = 1.0e-3 , a_x = 1.0 , a_dxdt = 1.0;
    controlled_stepper_type controlled_stepper( 
        default_error_checker< double , range_algebra , default_operations >( abs_err , rel_err , a_x , a_dxdt ) );

    integrate_adaptive( controlled_stepper , bar , x , std::log(10.0) , std::log(0.12), -0.1, push_back_state_and_time( x_vec, cutoffs) );


    std::cout << x_vec.size() << '\n';
    std::cout << cutoffs.size() << '\n';

    std::vector<double> zt(cutoffs.size());
    for(std::size_t i=0; i<cutoffs.size(); i++){
        zt[i] = x_vec[i][0];
    }
    
    seValuesReal sev;
    sev.m_val = cutoffs;
    sev.m_SEPos = zt;

    H5::H5File file(filenameTang, H5F_ACC_TRUNC );
    writeTang writefile(scale, bv, sev);
    writefile.writeMainResults(file);

    std::cout << "Results written to: \n " << filenameTang << '\n';

    double tseconds{ t.elapsed() };
    int minutes = (int) tseconds/60.0;
    int seconds = (int) tseconds - 60.0*minutes;
    
    std::cout << "Time taken: " << minutes << " minutes and " << seconds << " seconds\n";
    std::cout << "Program executed!\n";

    return 0;
}
