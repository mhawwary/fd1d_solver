#pragma once

#include "global_var.h"
#include "general_tools.h"

void intermediate_solution_dump(const int& IterNo, const double& time_);
void error_dumping(double& L1_norm_,double& L2_norm_);

void initial_solution_dumping(void);
void final_solution_dump(void);

void BinaryDataReading();
void BinaryDataWriting(const int& IterNo_);
void copy_sol_(double* Qn, double* Qn_temp1_);


//================================================================


void initial_solution_dumping(){

    register int i;

    FILE* sol_out=fopen("./output/u_initial.dat","w");

    for(i=Nghost_l;i<Nfaces+Nghost_l;i++)
    {
        fprintf(sol_out, "%2.10e\n", Qn[i]);
    }

    fclose(sol_out);

    FILE* coord_out=fopen("./output/x.dat","w");

    for(i=0;i<Nfaces;i++)
    {
        fprintf(coord_out, "%2.10e\n", X[i]);
    }

    fclose(coord_out);

    FILE* sol_out1=fopen("./output/initial_sol.dat","w");

    for(i=0;i<Nfaces;i++)
    {
        fprintf(sol_out1, "%2.10e %2.10e\n", X[i], Qn[i+Nghost_l]);
    }

    fclose(sol_out1);


    return;
}

void final_solution_dump(){

    register int i;

    char *fname=nullptr;
    fname = new char[100];

    //sprintf(fname,"./output/u_final_t%1.3f.dat",max_time);

    sprintf(fname,"./output/u_final.dat");

    FILE* sol_out=fopen(fname,"w");

    for(i=Nghost_l; i<Nfaces+Nghost_l; i++)
    {
        fprintf(sol_out, "%2.10e\n", Qn[i]);
    }

    fclose(sol_out);

    emptyarray(fname);

    fname = new char[100];

    sprintf(fname,"./output/u_exact.dat");

    FILE* sol_out1=fopen(fname,"w");

    for(i=0; i<Nfaces; i++)
    {
        fprintf(sol_out1, "%2.10e\n", Qex[i]);
    }

    fclose(sol_out1);

    emptyarray(fname);


    return;
}

void intermediate_solution_dump(const int& IterNo, const double& time_){


    return;
}

void error_dumping(double &L1_norm_, double &L2_norm_){

    char *error_fname=nullptr;
    error_fname= new char[100];

    sprintf(error_fname,"./output/error_CFL%1.3f_t%1.3f_N%d.dat",CFL,max_time,Nelem);

    FILE* error_out=fopen(error_fname,"at+");

    fprintf(error_out, "%d %e %e %e %e\n", Nelem, dt , CFL, L1_norm_, L2_norm_);

    fclose(error_out);

    emptyarray(error_fname);

    return;
}

void BinaryDataReading(){

    char *cfname=nullptr;  cfname=new char[200];

    sprintf(cfname,"./output/BINARY/u_binary_iter(%d).dat",restart_iter_);

    FILE*  bsol_read=fopen(cfname,"rb");

    fread(&Nelem,sizeof(int),1,bsol_read);
    fread(&scheme_order,sizeof(int),1,bsol_read);
    fread(&RK_order,sizeof(int),1,bsol_read);
    fread(&CFL,sizeof(double),1,bsol_read);

    fread(Qinit, sizeof(double),Nfaces,bsol_read);

    dx=(xf-x0)/Nelem;
    dt=dx*CFL/a_wave;

    t_init = dt * restart_iter_;

    fclose(bsol_read);

    emptyarray(cfname);

    return;
}

void BinaryDataWriting(const int& IterNo_){

    char *cfname=nullptr;  cfname=new char[200];

    sprintf(cfname,"./output/BINARY/u_binary_iter(%d).dat",IterNo_);


    FILE*  bsol_write=fopen(cfname,"wb");

    fwrite(&Nelem,sizeof(int),1,bsol_write);
    fwrite(&scheme_order,sizeof(int),1,bsol_write);
    fwrite(&RK_order,sizeof(int),1,bsol_write);
    fwrite(&CFL,sizeof(double),1,bsol_write);

    double* Qn_temp1_=nullptr;

    Qn_temp1_ = new double[Nfaces];

    copy_sol_(Qn, Qn_temp1_);

    fwrite(Qn_temp1_, sizeof(double),Nfaces,bsol_write);

    fclose(bsol_write);

    emptyarray(cfname);
    emptyarray(Qn_temp1_);

    return;
}

void copy_sol_(double* Qn, double* Qn_temp1_){

    register int i;

    for(i=0; i<Nfaces; i++)  Qn_temp1_[i] = Qn[i+Nghost_l];


    return;
}
