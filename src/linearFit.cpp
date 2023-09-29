#include "linearFit.h"

linearFit::linearFit(std::vector<double> x, std::vector<double> y): m_x(x), m_y(y) {}

// slope without constant term. i.e. y=f[0]=0
double linearFit::slope_noConstant(){
    const size_t n = m_x.size();

    double c1, cov11, chisq;

    gsl_fit_mul(&m_x[0],1,&m_y[0],1,n,&c1, &cov11, &chisq);
    return c1;
}

double linearFit::slope_regular(){
    const size_t n = m_x.size();

    double c0, c1, cov00, cov01, cov11, chisq;

    gsl_fit_linear(&m_x[0],1,&m_y[0],1,n, &c0, &c1, &cov00, &cov01, &cov11, &chisq);

    return c1;
}

double linearFit::intercept_regular(){
    const size_t n = m_x.size();

    double c0, c1, cov00, cov01, cov11, chisq;

    gsl_fit_linear(&m_x[0],1,&m_y[0],1,n, &c0, &c1, &cov00, &cov01, &cov11, &chisq);

    return c0;
}

