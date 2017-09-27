#ifndef PADEFILTER_H
#define PADEFILTER_H

#include"general_tools.h"
#include"global_var.h"
#include"solver_tools.h"

class PadeFilter{

public:

//  Construction Functions :
   PadeFilter(){}
   ~PadeFilter();

// member functions:

public:
   void setup_filter(const int& filter_OA_, const int& Nnodes_
                     , const double& alpha_filt_, string& bound_type_);
   void filtered_sol(double **qn_temp_);

private:
   void Reset_filter();
   void setup_filter_coeff();
   void update_unfiltered_sol(double** qq_t);
   void compute_RHS_vec(double *qn_temp_);
   void copy_Qfilt_to_sol(double **qn_temp_);

private:
   int filter_order=2;   // order of the filter
   int Nnodes=1;         // no. of nodes
   int n_linsys=1;
   int stencil_size=1;     // filter stencil size
   double alpha_f=0.0;     // alpha_filter
   double* alpha_f_vec=nullptr;
   double* d_vec=nullptr;
   double* RHS_vec = nullptr;
   double* coeff=nullptr;   // array of coeficients of the filter
   double* Q_filt=nullptr;  // filtered solution
   int Nghost_l=1;
   int Nghost_r=1;
   int Nnodes_tot=1;

   std::string Bound_type;

};

#endif
