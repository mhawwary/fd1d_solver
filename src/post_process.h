#pragma once

#include "global_var.h"
#include "general_tools.h"

void intermediate_solution_dump(const int& IterNo, const double& time_);
void error_dumping(double& L1_norm_,double& L2_norm_);

void initial_solution_dumping(void);
void initial_ghost_sol_dump();
void final_solution_dump(void);
//void refined_mesh_dumping(void);


//================================================================

void initial_ghost_sol_dump(){

    register int i;

    char *fname=nullptr; fname =new char[100];

    sprintf(fname, "./output/initial_ghost_sol_t%1.2e.dat",gtime);

    FILE* sol_out=fopen(fname,"w");

    for(i=0;i<Nfaces+Nghost_l+Nghost_r;i++)
    {
        fprintf(sol_out, "%2.10e\n", Qtemp[i]);
    }

    fclose(sol_out);

    return;
}

void initial_solution_dumping(){

    register int i;

    FILE* sol_out=fopen("./output/u_initial.dat","w");

    for(i=0;i<Nfaces;i++)
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
        fprintf(sol_out1, "%2.10e %2.10e\n", X[i], Qn[i]);
    }

    fclose(sol_out1);


    return;
}

void final_solution_dump(){

    register int i;

    char *fname=nullptr,*fname1=nullptr;
    fname = new char[100];
    fname1 = new char[100];

    //sprintf(fname,"./output/u_final_t%1.3f.dat",max_time);

    sprintf(fname,"./output/u_final.dat");

    FILE* sol_out=fopen(fname,"w");

    for(i=0; i<Nfaces; i++)
    {
        fprintf(sol_out, "%2.10e\n", Qn[i]);
    }

    fclose(sol_out);
    emptyarray(fname);

    sprintf(fname1,"./output/u_ghost.dat");

    FILE* sol_out1=fopen(fname1,"w");

    for(i=0; i<Nfaces+Nghost_l+Nghost_r; i++)
    {
        fprintf(sol_out1, "%2.10e\n", Qtemp[i]);
    }

    fclose(sol_out1);

    emptyarray(fname1);


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
