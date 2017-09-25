#ifndef FDSOLVERADVECDIFFUS_H
#define FDSOLVERADVECDIFFUS_H

#include "FDSolver.hpp"


class FDSolverAdvecDiffus:public FDSolver{

public:

//  Construction Functions :
   FDSolverAdvecDiffus(void);
   virtual ~FDSolverAdvecDiffus(void);

   virtual void setup_solver(GridData* meshdata_, SimData& simdata_);

   virtual void InitSol();
   virtual void UpdateResid(double **Resid_, double **Qn_);

   virtual double L1_error_nodal_sol();
   virtual double L2_error_nodal_sol();

   virtual void print_cont_vertex_sol();
   virtual void dump_errors(double &L1_error_, double &L2_error_);

   virtual void setup_coefficients();

   virtual void update_ghost_sol(double **Qn_);

   virtual double eval_init_sol(const double& xx);

   virtual void CalcTimeStep();

   virtual void Reset_solver();
   virtual void ComputeExactSolShift();

   virtual void Compute_exact_sol();
   virtual void Compute_exact_sol_for_plot();
   virtual void dump_timeaccurate_sol();

protected:
   double evaluate_inviscid_flux(const double& qn_);
   double eval_init_u_decay_burger_turb(const double& xx);
   void compute_RHS_f1_implicit(const double& hh_, double** qn_
                                , double*& RHS_temp_);
   void compute_RHS_f2_implicit(const double& hh_, double** qn_
                                , double*& RHS_temp_);

protected:
   // for explicit FD:
   double *FD_coeff_2nd=nullptr;  // coefficients for 2nd derivative
   int *stencil_index_2nd=nullptr; // stencil for 2nd derivative

   // for implicit FD:
   double *alpha_vec_f1_ =nullptr;  // vector of alpha_f1_, f'
   double *b_vec_ =nullptr;      // vector of 1 ones on each A matrix diagonal
   double *alpha_vec_f2_ =nullptr;  // vector of alpha_f2_, f"

   double *RHS_f1_=nullptr;
   double *RHS_f2_=nullptr;

   int n_linsys=0;

};

#endif
