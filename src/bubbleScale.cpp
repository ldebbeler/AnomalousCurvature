#include "bubbleScale.h"

bubbleScaleReal::bubbleScaleReal(double T): scaleEqns(T), m_spline(realSpline()) {}

double bubbleScaleReal::asymptoticLeft(double x) {
    return realAsymptotic(x);
}

double bubbleScaleReal::asymptoticRight(double x) {
    return realAsymptotic(x);
}

double bubbleScaleReal::evaluate(double x){
    if(x<-Lambda){
        return asymptoticLeft(x);
    }
    else if(x>Lambda){
        return asymptoticRight(x);
    }
    else{
        return m_spline.evaluate(x);
    }
}

bubbleScaleImag::bubbleScaleImag(double T): scaleEqns(T), m_spline(imagSpline()) {}

double bubbleScaleImag::asymptoticLeft(double x) {
    return imagAsymptotic(x);
}

double bubbleScaleImag::asymptoticRight(double x) {
    return imagAsymptotic(x);
}

double bubbleScaleImag::evaluate(double x){
    if(x<-Lambda){
        return asymptoticLeft(x);
    }
    else if(x>Lambda){
        return asymptoticRight(x);
    }
    else{
        return m_spline.evaluate(x);
    }
}


