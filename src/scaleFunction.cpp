#include "scaleFunction.h"
#include "constants.h"

scaleFunction::scaleFunction(interpolater1d spline): m_spline(spline) {}

double scaleFunction::evaluateSpline(double x){
    return evaluate(x);
}

double scaleFunction::evaluate(double x){
    if(x<-Lambda){
        return asymptoticLeft(x);
    }
    else if(x>Lambda){
        return asymptoticRight(x);
    }
    else{
        return evaluateSpline(x);
    }
}


