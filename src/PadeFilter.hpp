#ifndef PADEFILTER_H
#define PADEFILTER_H

#include"general_tools.h"
#include"global_var.h"
#include"solver_tools.h"
#include"filter.hpp"

class PadeFilter:public Filter{

public:

//  Construction Functions :
   PadeFilter(){}
   virtual ~PadeFilter();

// member functions:
   virtual void setup_filter(std::string& in_filter_type_,const int& filter_OA_
                             ,const int& Nnodes_,const double& alpha_filt_
                             ,string& bound_type_);
   virtual void filtered_sol(double **qn_temp_);

   virtual void Reset_filter();
   virtual void setup_filter_coeff();
   virtual void update_unfiltered_sol(double** qq_t);
   virtual void copy_Qfilt_to_sol(double **qn_temp_);

protected:
   void compute_RHS_vec(double *qn_temp_);

private:
   int n_linsys=1;
   double* alpha_f_vec=nullptr;
   double* d_vec=nullptr;
   double* RHS_vec = nullptr;
};

#endif
