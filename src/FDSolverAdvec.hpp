#ifndef FDSOLVERADVEC_H
#define FDSOLVERADVEC_H

#include "FDSolver.hpp"
#include "PadeFilter.hpp"


class FDSolverAdvec:public FDSolver{

public:

//  Construction Functions :
   FDSolverAdvec(void);
   virtual ~FDSolverAdvec(void);

   virtual void setup_solver(GridData* meshdata_, SimData& simdata_);

   virtual void InitSol();
   virtual void UpdateResid(double **Resid_, double **Qn_);

   virtual double L1_error_nodal_sol();
   virtual double L2_error_nodal_sol();
   virtual double dissipation_error();
   //virtual double L1_error_time_nodal_sol();
   //virtual double L2_error_time_nodal_sol();

   virtual void print_cont_vertex_sol();
   virtual void dump_errors(double &L1_error_, double &L2_error_);
   //virtual void dump_errors_vs_time(double &L1_error_, double &L2_error_);

   virtual void setup_coefficients();

   virtual void update_ghost_sol(double **Qn_);

   virtual double eval_init_sol(const double& xx);

   virtual void CalcTimeStep();

   virtual void Reset_solver();
   virtual void ComputeExactSolShift();

   virtual void Compute_exact_sol();
   virtual void Compute_exact_sol_for_plot();
   virtual void dump_timeaccurate_sol();
   virtual void dump_timeaccurate_errors();
   virtual void filter_solution(double **qn_);

protected:
   double evaluate_inviscid_flux(const double& qn_);
   void compute_RHS_f1_implicit(const double& hh_, double** qn_
                                , double*& RHS_temp_);
   void Compute_TimeAccurate_exact_sol();

protected:
   // for implicit FD:
   double *alpha_vec_f1_ =nullptr;  // vector of alpha_f1_, f'
   double *b_vec_ =nullptr;      // vector of 1 ones on each A matrix diagonal

   double *RHS_f1_=nullptr;

   PadeFilter *filter=nullptr;
   double filter_alpha_=0.0;
//   double *Qn_filt=nullptr; // one D array of nodal solutions

   int n_linsys=0;

   // for time accurate error calculations
   //double **qq_exact_time=nullptr;

};

#endif
