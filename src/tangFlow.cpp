#include "tangFlow.h"
#include "integrate2d.hpp"
#include "constants.h"

tangFlow::tangFlow(bubble foo): m_bubble(foo) {}

std::complex<double> tangFlow::fraction(double lambda, double qr, double qt, int num, int denom){
    return std::pow(qt,num)/std::pow(I*lambda+vF*qr-b*qt*qt*qt*qt,denom);
}

std::complex<double> tangFlow::fermionIntegrand(double lambda, double qr, double qt){
    return 256*b*b*b*fraction(lambda, qr, qt, 12, 5) + 288*b*b*fraction(lambda, qr, qt, 12, 5)
         + 68*b*fraction(lambda, qr, qt, 4, 3) + fraction(lambda,qr, qt,0,2);
}

std::complex<double> tangFlow::fullIntegrand(double lambda, double qr, double qt){
    return 1.0/(2.0*M_PI)*(fermionIntegrand(lambda, qr, qt) + fermionIntegrand(-lambda,qr,qt))*m_bubble.evaluateCombo(0.0,qr,qt);
}

std::complex<double> tangFlow::rightFlow(double lambda){
    auto integ = [&] (double2 Q) -> std::complex<double>{
        double Qr = Q.x;
        double Qt = Q.y;
        return std::exp(Qr)*std::exp(Qt)*(fullIntegrand(lambda, std::exp(Qr), std::exp(Qt))
                                         +fullIntegrand(lambda, std::exp(Qr),-std::exp(Qt))
                                         +fullIntegrand(lambda,-std::exp(Qr), std::exp(Qt))
                                         +fullIntegrand(lambda,-std::exp(Qr),-std::exp(Qt)));
    };
    return integrate2Dmeasure2pi(std::log(IR),std::log(UV),std::log(IR),std::log(UV)).integrate(integ,prec2d,steps2d);
}

void tangFlow::operator() (const std::vector<double>& x, std::vector<double>& dxdl, const double l){
    dxdl[0] = std::exp(l)*rightFlow(std::exp(l)).real();
}
