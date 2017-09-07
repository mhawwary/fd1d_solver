#ifndef FDSOLVER_H
#define FDSOLVER_H

#include"GridData.h"
#include"SimData.hpp"
#include"general_tools.h"
#include"global_var.h"

//struct SimData;
//struct GridData;

class FDSolver{

public:

//  Construction Functions :
   FDSolver(){}
   virtual ~FDSolver(){}

   virtual void setup_solver(GridData* meshdata_, SimData& simdata_)=0;

   virtual void InitSol()=0;
   virtual void UpdateResid(double **Resid_, double **Qn_)=0;
   //virtual void UpdateSolution(double **Qn_)=0;

   virtual double L1_error_nodal_sol()=0;
   virtual double L2_error_nodal_sol()=0;

   void UpdatePhyTime(const double& dt_){

       phy_time += dt_;

       return;
   }

   void SetPhyTime(const double &time_){

       phy_time=time_;

       return;
   }

   double GetTimeStep(){

       return time_step;
   }

   double GetLastTimeStep(){

       return last_time_step;
   }

   double GetCFL(){

       return CFL;
   }

   int GetNdof(){

       return Ndof;
   }

   double** GetNumSolution(){

       return Qn;
   }

   double** GetExactSolution(){

       return Q_exact_pp;
   }

   double GetPhyTime(){

       return phy_time;
   }

   int GetNfaces(){
       return Nfaces;
   }

   int GetNghost_l(){
       return Nghost_l;
   }

   virtual void print_cont_vertex_sol()=0;
   virtual void dump_errors(double &L1_error_, double &L2_error_)=0;

protected:

   virtual void setup_stencil()=0;

   virtual void update_ghost_sol(double **Qn_)=0;

   virtual double eval_init_sol(const double& xx)=0;

   virtual void CalcTimeStep()=0;

   virtual void Reset_solver()=0;
   virtual void ComputeExactSolShift()=0;

   virtual void Compute_exact_sol()=0;
   virtual void Compute_exact_sol_for_plot()=0;

protected:

   GridData *grid_=nullptr;
   SimData *simdata_=nullptr;

   int scheme_order_=1;
   int *stencil_index=nullptr;
   double *FD_coeff=nullptr;

   int Ndof = 1;

   double **Qn=nullptr;      // Nfaces_tot * Ndof long

   double **Qn_true=nullptr; // Nfaces * Ndof long

   double **Q_init=nullptr;  // Nfaces * Ndof long

   double **Q_exact_pp=nullptr; //Nppoints * Ndof long

   double **Q_exact=nullptr;

   double phy_time=0.0;
   double time_step=1e-5;
   double last_time_step=1e-5;
   double CFL=1.0;

   double exact_sol_shift=0.;
   double wave_length_=0.;
   double T_period=1;

   int Nfaces=1;
   int Nghost_l=1;
   int Nghost_r=1; // ghost points at the left and right interfaces
   int Nfaces_tot=1;

};

#endif
