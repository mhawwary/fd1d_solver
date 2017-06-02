#include"GridData.h"

void GridData::set_grid_param(const SimData& simdata_){

    Nelem = simdata_.Nelem_;

    x0 = simdata_.x0_;
    xf = simdata_.xf_;

    uniform = simdata_.uniform_;

    refine_level=simdata_.refine_level_;

    Nfaces = Nelem+1;

    X = new double[Nfaces];

    h_j = new double [Nelem];

    N_exact_ppts= simdata_.N_exact_ppts;

    x_exact_ppts = new double[N_exact_ppts];

    return;
}

void GridData::generate_grid(){

    register int i;
    
    if(uniform==1) dx = (xf-x0)/Nelem;
    else if(uniform==0)
        _notImplemented("non-unifrom grid treatment");

    Idx = 1.0/dx;

    for(i=0; i<Nfaces; i++)
        X[i]   = dx * (i)  + x0 ;  // node 0, element i

    for (i=0; i<Nelem; i++)
        h_j[i]= X[i+1]-X[i];

    // New sampling for plotting a smooth exact solution

    double dxx_ = (xf - x0) / (N_exact_ppts-1);

    for(i=0; i<N_exact_ppts; i++){

        x_exact_ppts[i]   = dxx_ * (i)  + x0 ;
    }

    return;
}
