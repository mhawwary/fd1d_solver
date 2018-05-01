#include "explicitfilter.hpp"

// Constructor/Destructor/ Setup functions:

ExplicitFilter::~ExplicitFilter()
{
    Reset_filter();
}

void ExplicitFilter::setup_filter(std::string& in_filter_type_,const int& stencil_size_
                                  ,const int& Nnodes_,const double& alpha_filt_
                                  ,string& bound_type_){

    Bound_type = bound_type_;  // Current implementation allows only for periodic boundary conditions
    stencil_size = stencil_size_;
    filter_type= in_filter_type_;
    Nnodes = Nnodes_;   // no. of true mesh nodes
    alpha_f = alpha_filt_;

    stencil_oneside_size = (stencil_size-1)/2+1;
    Nghost_l = stencil_oneside_size-1;
    Nghost_r = stencil_oneside_size-1;
    Nnodes_tot = Nghost_l + Nnodes + Nghost_r;

    printf("\nfilter_order:%s\nfilter_stencil_size:%d\nstencil_halfsize:%d\nalpha_filter:%1.3f\n"
           ,filter_type.c_str(), stencil_size, stencil_oneside_size, alpha_f);

    stencil_index = new int[stencil_oneside_size];
    coeff = new double[stencil_oneside_size];
    Q_filt = new double[Nnodes_tot];

    setup_filter_coeff();

    return;
}

void ExplicitFilter::Reset_filter(){

    emptyarray(stencil_index);
    emptyarray(coeff);
    emptyarray(Q_filt);

    return;
}

// Member functions
//-------------------------------------------
void ExplicitFilter::setup_filter_coeff(){

    register int i;

    for(i=0; i<stencil_oneside_size; i++)  stencil_index[i] = i;

    switch (stencil_size){
    case 9:  // SF9
        if(filter_type=="BogeyBailly"){
            coeff[0] =  0.243527493120;
            coeff[1] = -0.204788880640;
            coeff[2] =  0.120007591680;
            coeff[3] = -0.045211119360;
            coeff[4] =  0.008228661760;
        }else if(filter_type=="standard") { // Standard explicit
            coeff[0] =  35./128.;
            coeff[1] = -7./32.;
            coeff[2] =  7./64.;
            coeff[3] = -1./32;
            coeff[4] =  1./256.;
        }
        break;

    case 11:  // SF11
        if(filter_type=="BogeyBailly"){
            coeff[0] =  0.215044884112;
            coeff[1] = -0.187772883589;
            coeff[2] =  0.123755948787;
            coeff[3] = -0.059227575576;
            coeff[4] =  0.018721609157;
            coeff[5] = -0.002999540835;
        }else if(filter_type=="standard") { // Standard explicit
            coeff[0] =  63./256.;
            coeff[1] = -105./512.;
            coeff[2] =  15./128.;
            coeff[3] = -45./1024.;
            coeff[4] =  5./512.;
            coeff[5] = -1./1024.;
        }
        break;

    case 13:  // SF13
        if(filter_type=="BogeyBailly"){
            coeff[0] =  0.190899511506;
            coeff[1] = -0.171503832236;
            coeff[2] =  0.123632891797;
            coeff[3] = -0.069975429105;
            coeff[4] =  0.029662754736;
            coeff[5] = -0.008520738659;
            coeff[6] =  0.001254597714;
        }else if(filter_type=="standard") { // Standard explicit
            coeff[0] =  231./1024.;
            coeff[1] = -99./512.;
            coeff[2] =  495./4096.;
            coeff[3] = -55./1024.;
            coeff[4] =  33./2048.;
            coeff[5] = -3./1024.;
            coeff[6] =  1./4096.;
        }
        break;

        break;
    default:
        break;
    }

    return;
}

void ExplicitFilter::update_unfiltered_sol(double** qq_t){

    register int i;

    for(i=0; i< Nnodes; i++)
        Q_filt[i+Nghost_l] = qq_t[i][0];

    // updating the left boundary:
    for(i=0; i<Nghost_l; i++)
        Q_filt[i] = Q_filt[Nnodes-1+i];

    //updating the right boundary:
    for(i=0; i<Nghost_r; i++)
        Q_filt[i+Nghost_l+Nnodes] = Q_filt[i+Nghost_l+1];

    return;
}

void ExplicitFilter::filtered_sol(double** qn_temp_){

    update_unfiltered_sol(qn_temp_);

    //perform filtering and excludes boundary points

    register int i=0; int j,s1;
    double II=0.0;

    for(i=1; i<Nnodes-1; i++){ // excludes boundary points
        II = Q_filt[i+Nghost_l+0]*coeff[0];  // u0
        for(j=1; j<stencil_oneside_size; j++){
            s1 = stencil_index[j];
            II += (Q_filt[i+Nghost_l+s1]
                    + Q_filt[i+Nghost_l-s1])*coeff[j] ;
        }
        qn_temp_[i][0] = Q_filt[i+Nghost_l+0]-alpha_f*II;
        II=0.0;
    }

    //copy_Qfilt_to_sol(qn_temp_);

    return ;
}

void ExplicitFilter::copy_Qfilt_to_sol(double **qn_temp_){

    register int i;

    for(i=1; i<Nnodes-1; i++)
        qn_temp_[i][0] = Q_filt[i+Nghost_l];

    return;
}
