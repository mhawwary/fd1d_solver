#ifndef FILTER_HPP
#define FILTER_HPP

#include"general_tools.h"
#include"global_var.h"
#include"solver_tools.h"

class Filter
{
public:

//  Construction Functions :
    Filter(){}
    void setFilterType(std::string in_filter_type){
        filter_type=in_filter_type;
    }

    virtual ~Filter(){}

// member functions:
    virtual void setup_filter(std::string& filter_type_, const int& filter_driver_
                      ,const int& Nnodes_,const double& alpha_filt_
                      , string& bound_type_)=0;
    //filter_driver is either filter_order for Pade or stencil_size for explicit
    virtual void filtered_sol(double **qn_temp_)=0;

protected:
    virtual void Reset_filter()=0;
    virtual void setup_filter_coeff()=0;
    virtual void update_unfiltered_sol(double** qq_t)=0;
    virtual void copy_Qfilt_to_sol(double **qn_temp_)=0;

protected:
   std::string filter_type; // Pade or Explicit == standard or BogeyBailly2004/BogeyBailly2006
   int filter_order=2;   // order of the filter
   int Nnodes=1;         // no. of nodes
   int stencil_size=1;     // filter stencil size
   double alpha_f=0.0;     // alpha_filter
   double* coeff=nullptr;   // array of coeficients of the filter
   double* Q_filt=nullptr;  // filtered solution
   int Nghost_l=1;
   int Nghost_r=1;
   int Nnodes_tot=1;

   std::string Bound_type;
};

#endif // FILTER_HPP
