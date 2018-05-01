#ifndef EXPLICITFILTER_HPP
#define EXPLICITFILTER_HPP

#include"general_tools.h"
#include"global_var.h"
#include"solver_tools.h"
#include"filter.hpp"

class ExplicitFilter:public Filter{

public:

//  Construction Functions :
   ExplicitFilter(){}
   virtual ~ExplicitFilter();

// member functions:

public:
   virtual void setup_filter(std::string& in_filter_type_,const int& stencil_size_
                             ,const int& Nnodes_,const double& alpha_filt_
                             ,string& bound_type_);
   virtual void filtered_sol(double **qn_temp_);

   virtual void Reset_filter();
   virtual void setup_filter_coeff();
   virtual void update_unfiltered_sol(double** qq_t);
   virtual void copy_Qfilt_to_sol(double **qn_temp_);

private:
   int stencil_oneside_size=1; // half stencil width
   int* stencil_index=nullptr;

};

#endif // EXPLICITFILTER_HPP
