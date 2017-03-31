
void setup_stencil();

void ComputeOneStep();

void compute_residual();

void update_sol();


//===========================================================

void setup_stencil(){

    if(scheme_order==1){      // first order upwind scheme

        stencil_index[0] =  0;
        stencil_index[1] = -1;           //[j,j-1];

        FD_coeff [0] =  1;
        FD_coeff [1] = -1;

    } else if (scheme_order==2) { // 2nd order central scheme

        stencil_index[0] =  1;
        stencil_index[1] =  0;
        stencil_index[2] = -1;           //[j+1,j,j-1];

        FD_coeff [0] =  0.5;
        FD_coeff [1] =  0.0;
        FD_coeff [2] = -0.5;

    }else if (scheme_order==3){

        // for 3rd order fully upwind:
        //stencil_index[0] =  0;
        //stencil_index[1] = -1;
        //stencil_index[2] = -2;
        //stencil_index[3] = -3;  // [j,j-1,j-2,j-3];

        // for 3rd order biased upwind:
        stencil_index[0] =  1;
        stencil_index[1] =  0;
        stencil_index[2] = -1;
        stencil_index[3] = -2;  //[j+1,j,j-1,j-2];

        FD_coeff [0] =  1.0/3.0;
        FD_coeff [1] =  0.5;
        FD_coeff [2] = -1.0;
        FD_coeff [3] =  1.0/6.0;

    }else if (scheme_order==4){ // 4th order central scheme

        stencil_index[0] =  2;
        stencil_index[1] =  1;
        stencil_index[2] =  0;
        stencil_index[3] = -1;
        stencil_index[4] = -2;  //[j+2,j+1,j,j-1,j-2];

        FD_coeff [0] =  -1.0/12.0;
        FD_coeff [1] =   2.0/3.0;
        FD_coeff [2] =   0.0;
        FD_coeff [3] =  -2.0/3.0;
        FD_coeff [4] =   1.0/12.0;
    }

    return;
}

void compute_residual(){




    return;
}

void compute_residual(){

    register int i,j;

    int s=0;

    double temp=0.0;

    for(i=1; i<Nfaces-1; i++){

        temp=0.0;

        for(j=0; j<scheme_order+1; j++){

            s = stencil_index[j];
            temp += Qn[i+s] * FD_coeff[j];
        }

        Resid[i] = CFL * temp;
    }


    return;

}


void ComputeOneStep(){

    compute_residual();

    time_integrate();

    update_sol();

    return;
}
