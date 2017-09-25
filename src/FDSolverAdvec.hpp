#ifndef FDSOLVERADVEC_H
#define FDSOLVERADVEC_H

#include "FDSolver.hpp"


class FDSolverAdvec:public FDSolver{

public:

//  Construction Functions :
   FDSolverAdvec(void);
   virtual ~FDSolverAdvec(void);

   virtual void setup_solver(GridData* meshdata_, SimData& simdata_);

   virtual void InitSol();
   virtual void UpdateResid(double **Resid_, double **Qn_);
   //virtual void UpdateSolution(double **Qn_);

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
   virtual void   dump_timeaccurate_sol();

};

#endif
