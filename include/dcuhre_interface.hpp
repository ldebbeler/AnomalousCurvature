# pragma once

extern "C" {

    typedef void (*integrand_dcuhre)(double*,int,double*,int, void*);

    // Fortran subroutine
    extern void dcuhre_interface_2D(int funcSize, double x1, double x2, double y1, double y2, 
        double acc_rel, double acc_abs, double *Res, integrand_dcuhre, void*, int *err);

}
