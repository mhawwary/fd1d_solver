#include"PadeFilter.hpp"

// Constructor/Destructor/ Setup functions:
//------------------------------------------------
PadeFilter::~PadeFilter(){
    Reset_filter();
}

void PadeFilter::setup_filter(const int& filter_OA_, const int& Nnodes_
                              ,const double& alpha_filt_, string& bound_type_){

    Bound_type = bound_type_;  // Current implementation allows only for periodic boundary conditions
    filter_order = filter_OA_;
    Nnodes = Nnodes_;   // no. of true mesh nodes
    n_linsys = Nnodes-2; // exculding the boundary points ( left and right )
    alpha_f = alpha_filt_;

    Nghost_l = (filter_order)/2-1;
    Nghost_r = (filter_order)/2-1;
    Nnodes_tot = Nghost_l + Nnodes + Nghost_r;

    stencil_size = (filter_order)/2 + 1;

    printf("\nfilter_order:%d\nfilter_stencil_size:%d\nalpha_filter:%1.3f\n"
           ,filter_order, stencil_size, alpha_f);

    coeff = new double[stencil_size];
    Q_filt = new double[Nnodes_tot];

    alpha_f_vec = new double[n_linsys];
    d_vec = new double[n_linsys];
    RHS_vec = new double[n_linsys];

    setup_filter_coeff();

    return;
}

void PadeFilter::Reset_filter(){

    emptyarray(alpha_f_vec);
    emptyarray(d_vec);
    emptyarray(RHS_vec);
    emptyarray(coeff);
    emptyarray(Q_filt);

    return;
}

// Member functions
//-------------------------------------------
void PadeFilter::setup_filter_coeff(){

    register int i;

    for(i=0; i<n_linsys; i++){
        alpha_f_vec[i] = alpha_f;
        d_vec[i] = 1.0;
        RHS_vec[i]=0.0;
    }

    switch (filter_order) {
    case 2:
        coeff[0] = 0.5 + alpha_f;
        coeff[1] = 0.5 * coeff[0];
        break;

    case 4:
        coeff[0] = (5./8.) + (3.*alpha_f/4.);
        coeff[1] = 0.5 * (0.5 + alpha_f);
        coeff[2] = 0.5 * ((-1./8.) + (alpha_f/4.));
        break;

    case 6:
        coeff[0] = (11./16.) + (5.*alpha_f/8.);
        coeff[1] = 0.5 * ((15./32.) + (17.*alpha_f/16.));
        coeff[2] = 0.5 * ((-3./16.) + (3.*alpha_f/8.));
        coeff[3] = 0.5 * ((1./32.) - (alpha_f/16.));
        break;

    case 8:
        coeff[0] = (93.+(70.*alpha_f))/128.;
        coeff[1] = 0.5 * (7.+(18.*alpha_f))/16.;
        coeff[2] = 0.5 * (-7.+(14.*alpha_f))/32.;
        coeff[3] = 0.5 * (0.0625 - (0.125*alpha_f));
        coeff[4] = 0.5 * (-1.+(2.*alpha_f))/128.;
        break;

    case 10:
        coeff[0] = (193.+(126.*alpha_f))/256.;
        coeff[1] = 0.5 * (105.+(302.*alpha_f))/256.;
        coeff[2] = 0.5 * (-15.+(30.*alpha_f))/64.;
        coeff[3] = 0.5 * (45.+(-90.*alpha_f))/512.;
        coeff[4] = 0.5 * (-5.+(10.*alpha_f))/256.;
        coeff[5] = 0.5 * (1.+(-2.*alpha_f))/512.;
        break;
    default:
        break;
    }

    return;
}

void PadeFilter::update_unfiltered_sol(double** qq_t){

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

void PadeFilter::compute_RHS_vec(double *qn_temp_){

    register int i;
    int j;

    for(i=0; i<n_linsys; i++){

        RHS_vec[i] = coeff[0] * qn_temp_[i+1];  // starting from q1

        for(j=1; j<stencil_size; j++)
            RHS_vec[i] += coeff[j] * (qn_temp_[i+j+1]+ qn_temp_[i-j+1]);
    }

    RHS_vec[0] += (- alpha_f * qn_temp_[0]);
    RHS_vec[n_linsys-1] += (- alpha_f * qn_temp_[n_linsys+1]);

    return;
}

void PadeFilter::filtered_sol(double** qn_temp_){

    update_unfiltered_sol(qn_temp_);
    compute_RHS_vec(&Q_filt[Nghost_l]);

    tridiag_solve_mh(n_linsys,alpha_f_vec,d_vec
                            ,alpha_f_vec,RHS_vec,&Q_filt[Nghost_l+1]);

    copy_Qfilt_to_sol(qn_temp_);

    return ;
}

void PadeFilter::copy_Qfilt_to_sol(double **qn_temp_){

    register int i;

    for(i=1; i<Nnodes-1; i++){
        qn_temp_[i][0] = Q_filt[i+Nghost_l];
    }

    return;
}




