#include "testBubble.h"
#include <iostream>

testBubble::testBubble(bubble foo): m_bubble(foo) {}

void testBubble::writeBubble(scalingValues& scale, bubbleValues& bv){
    scale.m_x = m_bubble.getxVector();
    scale.m_real= m_bubble.getReals();
    scale.m_imag= m_bubble.getImags();
    scale.m_T = m_bubble.getExp();

    std::vector<double> bubbleOmegaReal(201);
    std::vector<double> bubbleOmegaImag(201);
    std::vector<double> omegas(201);
    double omegqr{ 0.05 };
    double omegqt{ 0.4 };
    for(std::size_t i=0; i<omegas.size(); i++){
        double val{ -1.0 + 0.01*i };
        omegas[i] = val;
        bubbleOmegaReal[i] = m_bubble.evaluate(val,omegqr,omegqt).real();
        bubbleOmegaImag[i] = m_bubble.evaluate(omegas[i],omegqr,omegqt).imag();
    }

    std::vector<double> bubbleRadialReal(201);
    std::vector<double> bubbleRadialImag(201);
    std::vector<double> radials(201);
    double radOmega{ 1e-5 };
    double radqt{ 1e-5 };
    for(std::size_t i=0; i<radials.size(); i++){
        radials[i] = -1.0+0.01*i;
        bubbleRadialReal[i] = m_bubble.evaluate(radOmega,radials[i],radqt).real();
        bubbleRadialImag[i] = m_bubble.evaluate(radOmega,radials[i],radqt).imag();
    }

    std::vector<double> bubbleTangentialReal(100);
    std::vector<double> bubbleTangentialImag(100);
    std::vector<double> tangentials(100);
    double tangOmega{ 1e-5 };
    double tangqr{ 2e-5 };
    for(std::size_t i=0; i<tangentials.size(); i++){
        tangentials[i] = 0.01*(i+1);
        bubbleTangentialReal[i] = m_bubble.evaluate(tangOmega,tangqr,tangentials[i]).real();
        bubbleTangentialImag[i] = m_bubble.evaluate(tangOmega,tangqr,tangentials[i]).imag();
    }

    //bubbleValues bv;
    bv.m_omega = omegas;
    bv.m_omegReal = bubbleOmegaReal;
    bv.m_omegImag = bubbleOmegaImag;
    bv.m_omegqr = omegqr;
    bv.m_omegqt = omegqt;
    bv.m_qr = radials;
    bv.m_radReal = bubbleRadialReal;
    bv.m_radImag = bubbleRadialImag;
    bv.m_radOmega = radOmega;
    bv.m_radqt = radqt;
    bv.m_qt = tangentials;
    bv.m_tangReal = bubbleTangentialReal;
    bv.m_tangImag = bubbleTangentialImag;
    bv.m_tangOmega = tangOmega;
    bv.m_tangqr = tangqr;
}
 
