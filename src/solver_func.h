#pragma once

#include "global_var.h"
#include "general_tools.h"

void setup_stencil();

void ComputeOneStep();

void compute_residual();
void compute_residual1();

void update_ghost_sol();

void time_integrate();
void fwdEuler();
void SSPRK22();
void SSPRK33();
void SSPRK44();



//=======================================================================

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

void compute_residual1(){

    register int i,j;

    int s=0;

    double temp=0.0;

    if(scheme_order==1){

        temp += Qn[0] * FD_coeff[0];

        temp += Qn[Nfaces-2] * FD_coeff[1];

        Resid[0] = -temp/dx;

        for(i=1; i<Nfaces; i++){

            temp=0.0;

            for(j=0; j<scheme_order+1; j++){

                s = stencil_index[j];
                temp += Qn[i+s] * FD_coeff[j];
            }

            Resid[i] =  - temp / dx;
        }


    }else if (scheme_order==2){

        i=0; temp=0.0;

        for(j=0; j<scheme_order; j++){

            s = stencil_index[j];
            temp += Qn[i+s] * FD_coeff[j];
        }

        temp += Qn[Nfaces-2] * FD_coeff[2];

        Resid[i] = -temp/dx;

        for(i=1; i<Nfaces-1; i++){

            temp=0.0;

            for(j=0; j<scheme_order+1; j++){

                s = stencil_index[j];
                temp += Qn[i+s] * FD_coeff[j];
            }

            Resid[i] =  - temp / dx;
        }

        i=Nfaces-1;

        temp=0.0;

        temp += Qn[1] * FD_coeff[0];

        for(j=1; j<scheme_order+1; j++){

            s = stencil_index[j];
            temp += Qn[i+s] * FD_coeff[j];
        }

        Resid[i] = -temp/dx;



    }else if(scheme_order==3){

        i=0; temp=0.0;

        for(j=0; j<scheme_order-1; j++){

            s = stencil_index[j];
            temp += Qn[i+s] * FD_coeff[j];
        }

        temp += Qn[Nfaces-2] * FD_coeff[2];

        temp += Qn[Nfaces-3] * FD_coeff[3];

        Resid[i] = -temp/dx;

        for(i=1; i<Nfaces-1; i++){

            temp=0.0;

            for(j=0; j<scheme_order+1; j++){

                s = stencil_index[j];
                temp += Qn[i+s] * FD_coeff[j];
            }

            Resid[i] =  - temp / dx;
        }

        i=Nfaces-1;

        temp=0.0;

        temp += Qn[1] * FD_coeff[0];

        for(j=1; j<scheme_order+1; j++){

            s = stencil_index[j];
            temp += Qn[i+s] * FD_coeff[j];
        }

        Resid[i] = -temp/dx;

    }else if(scheme_order==4){

        i=0; temp=0.0;

        for(j=0; j<3; j++){

            s = stencil_index[j];
            temp += Qn[i+s] * FD_coeff[j];
        }

        temp += Qn[Nfaces-2] * FD_coeff[3];
        temp += Qn[Nfaces-3] * FD_coeff[4];

        Resid[i] = -temp/dx;

        for(i=1; i<Nfaces-1; i++){

            temp=0.0;

            for(j=0; j<scheme_order+1; j++){

                s = stencil_index[j];
                temp += Qn[i+s] * FD_coeff[j];
            }

            Resid[i] =  - temp / dx;
        }

        i=Nfaces-1;

        temp=0.0;

        temp += Qn[1] * FD_coeff[0];
        temp += Qn[2] * FD_coeff[1];

        for(j=2; j<scheme_order+1; j++){

            s = stencil_index[j];
            temp += Qn[i+s] * FD_coeff[j];
        }

        Resid[i] = -temp/dx;

    }


    return;
}

void compute_residual(){

    register int i,j;

    int s=0;

    double temp=0.0;

    for(i=Nghost_l; i<Netot-1; i++){

        temp=0.0;

        for(j=0; j<scheme_order+1; j++){

            s = stencil_index[j];
            temp += Qn[i+s] * FD_coeff[j];
        }

        Resid[i-Nghost_l] =  - temp / dx;

        //Resid[i-Nghost_l] = - CFL *(Qn[i]-Qn[i-1]);
    }

    if(scheme_order==4){

        temp=0.0;
        i = Netot-1;

        for(j=2; j<scheme_order+1; j++){

            s = stencil_index[j];
            temp += Qn[i+s] * FD_coeff[j];
        }

        temp += Qn[Nghost_l+2] *  FD_coeff[0];
        temp += Qn[Nghost_l+1] *  FD_coeff[1];

        Resid[i-Nghost_l] =  - temp / dx;

        temp=0.0;

        i=Netot-2;

        for(j=1; j<scheme_order+1; j++){

            s = stencil_index[j];
            temp += Qn[i+s] * FD_coeff[j];
        }

        temp += Qn[Nghost_l+1] *  FD_coeff[0];

        Resid[i-Nghost_l] =  - temp / dx;



    }else if(scheme_order==3){
        temp=0.0;
        i = Netot-1;

        for(j=1; j<scheme_order+1; j++){

            s = stencil_index[j];
            temp += Qn[i+s] * FD_coeff[j];
        }

        temp += Qn[Nghost_l+1] *  FD_coeff[0];

        Resid[i-Nghost_l] =  - temp / dx;


    }else if(scheme_order==2){

        temp=0.0;
        i = Netot-1;

        for(j=1; j<scheme_order+1; j++){

            s = stencil_index[j];
            temp += Qn[i+s] * FD_coeff[j];
        }

        temp += Qn[Nghost_l+1] *  FD_coeff[0];

        Resid[i-Nghost_l] =  - temp / dx;

    }else if(scheme_order==1){

        temp=0.0;
        i = Netot-1;

        for(j=0; j<scheme_order+1; j++){

            s = stencil_index[j];
            temp += Qn[i+s] * FD_coeff[j];
        }

        Resid[i-Nghost_l] =  - temp / dx;
    }

    return;

}

void update_ghost_sol(){

    register int i;

    for(i=0; i<Nghost_l; i++){

        Qn[i] = Qn[Netot-Nghost_l-1+i];
    }

    //Qn[0] = Qn[Netot-2];

//    for(i=0; i<Nghost_r; i++){

//        Qn[i] = Qn[Nghost_l+1+i];
//    }

    return;
}

void ComputeOneStep(){

    //update_ghost_sol();

    compute_residual1();

    time_integrate();

    return;
}


//=======================================================================


//=======================================================================

void time_integrate(){

    if(RK_order==1){

        fwdEuler();

    }else if(RK_order==2){

        SSPRK22();

    } else if(RK_order==3){

        SSPRK33();

    } else if(RK_order==4){

        SSPRK44();
    }

    return;
}

void fwdEuler(){

    register int i;

    for(i=0; i<Nfaces; i++){

        Qn[i] = Qn[i] + dt * Resid[i];
    }

    //Qn[Netot-1] = Qn[Nghost_l];

    return;
}

void SSPRK22(){

    register int i;

    double *q_temp=nullptr;
    q_temp = new double[Nfaces];

    for(i=0;i<Nfaces;i++)  q_temp[i]=Qn[i];

    // Step1:
    //-----------

    for(i=0; i<Nfaces; i++){

        Qn[i] = q_temp[i] + dt * Resid[i];
    }

    //Qn[Netot-1] = Qn[Nghost_l];

    //update_ghost_sol();

    compute_residual1();

    // Step2:
    //------------

    for(i=0; i<Nfaces; i++){

        Qn[i] = 0.5 * ( q_temp[i] +  Qn[i] + dt * Resid[i] );
    }

    //Qn[Netot-1] = Qn[Nghost_l];

    emptyarray(q_temp);

    return;
}

void SSPRK33(){

    register int i;

    double *q_temp=nullptr;
    q_temp = new double[Nfaces];

    for(i=0;i<Nfaces;i++)  q_temp[i]=Qn[i];

    // Step1:
    //-----------
    for(i=0; i<Nfaces; i++){

        Qn[i] = q_temp[i] + dt * Resid[i];
    }

    //Qn[Netot-1] = Qn[Nghost_l];

    //update_ghost_sol();

    compute_residual1();

    // Step2:
    //------------

    for(i=0; i<Nfaces; i++){

        Qn[i] = ( 0.75 *  q_temp[i] ) + 0.25 *( Qn[i] + dt * Resid[i] ) ;
    }

    //Qn[Netot-1] = Qn[Nghost_l];

    //update_ghost_sol();

    compute_residual1();

    //for(i=0;i<Netot;i++)  q_temp[i]=Qn[i];

    // Step3:
    //--------------

    for(i=0; i<Nfaces; i++){

        Qn[i] = ( (1.0/3.0) *  q_temp[i] )
                + (2.0/3.0) *( Qn[i] + dt * Resid[i] ) ;
    }

    //Qn[Netot-1] = Qn[Nghost_l];

    emptyarray(q_temp);

    return;
}

void SSPRK44(){

    return;
}


//=====================================================================
void ComputeError(double &L1_norm_, double &L2_norm_);
void Compute_exact_shifted_sol();
//=====================================================================

void Compute_exact_shifted_sol(){

    register int i;

    for(i=0; i<Nfaces; i++){

        Qex[i] = sin(2*PI*(X[i]-max_time*a_wave));
    }

    return;
}

void ComputeError(double &L1_norm_, double &L2_norm_){

    Compute_exact_shifted_sol();

    register int i;

    double l1=0.0,l2=0.0;

    for(i=0; i<Nfaces; i++){

        l1 += fabs(Qex[i]-Qn[Nghost_l+i]);
        l2 += pow(Qex[i]-Qn[Nghost_l+i],2);
    }

    L1_norm_ = l1 / Nfaces;
    L2_norm_ = l2 / Nfaces;

    return;
}
