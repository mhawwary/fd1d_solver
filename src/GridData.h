#ifndef GRIDDATA_H
#define GRIDDATA_H

#include "general_tools.h"
#include"SimData.hpp"
#include"../include/error.h"

//struct SimData;

struct GridData {

    double x0=0;
    double xf=1.0;
    double dx=1.0;
    double Idx=1.0;  // Idx = 1/dx;

    double *X=nullptr;    // mesh X coord
    double *h_j=nullptr;  // mesh width for each elements

    double *x_exact_ppts=nullptr;
    int N_exact_ppts=100;

    int Nelem=1;
    int Nfaces=2;

    int uniform=1;  // 0: for nonuniform mesh elements

    int refine_level=0; // 0: no refinement

    void set_grid_param(const SimData& simdata_);
    void generate_grid();

    ~GridData(){
        emptyarray(X);
        emptyarray(h_j);
        emptyarray(x_exact_ppts);
    }
};


#endif
