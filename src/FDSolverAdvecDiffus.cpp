#include "FDSolverAdvecDiffus.hpp"

// Constructor/Destructor/ Setup functions:
//------------------------------------------------
FDSolverAdvecDiffus::FDSolverAdvecDiffus(void){}

FDSolverAdvecDiffus::~FDSolverAdvecDiffus(void){

    Reset_solver();
}

void FDSolverAdvecDiffus::setup_solver(GridData* meshdata_, SimData& osimdata_){

    grid_ = meshdata_;
    simdata_ = &osimdata_;

    Ndof= 1;

    Nfaces = grid_->Nfaces;

    scheme_order_ = simdata_->scheme_order_;

    if(scheme_order_==1){
        Nghost_l = 1;
        Nghost_r=0;

    }else if(scheme_order_==2){
        Nghost_l=1;
        Nghost_r=1;

    }else if(scheme_order_==3){

        if(simdata_->upwind_biased_==1){
            Nghost_l=2;
            Nghost_r=1;
        }else if(simdata_->upwind_biased_==0){
            Nghost_l=3;
            Nghost_r=0;
        }

    }else if(scheme_order_==4){
        Nghost_l=2;
        Nghost_r=2;

    }else if(scheme_order_==6){
        Nghost_l=3;
        Nghost_r=3;
    }

    Nfaces_tot = Nghost_l + Nfaces + Nghost_r;

    Qn      =  new double* [Nfaces_tot];
    Q_init  =  new double* [Nfaces];
    Q_exact =  new double* [Nfaces];
    Q_exact_pp = new double*[grid_->N_exact_ppts];

    register int i;

    for(i=0; i<Nfaces_tot; i++)
        Qn[i]     = new double[Ndof];

    for(i=0; i<Nfaces; i++){
        Q_init[i] = new double[Ndof];
        Q_exact[i] = new double[Ndof];
    }

    for(i=0; i<grid_->N_exact_ppts; i++)
        Q_exact_pp[i] = new double [Ndof];

    setup_stencil();

    SetPhyTime(simdata_->t_init_);
    CalcTimeStep();
    ComputeExactSolShift();
    Compute_exact_sol();
    Compute_exact_sol_for_plot();

    // Screen Output of input and simulation parameters:
    cout <<"\n===============================================\n";
    cout << "CFL no.        : "<<CFL<<endl;
    cout << "time step, dt  : "<<time_step<<endl;
    cout << "last_time_step: "<<last_time_step<<endl;
    cout << "input Nperiods : "<<simdata_->Nperiods<<endl;
    cout << "new   Nperiods : "<<simdata_->t_end_/T_period<<endl;
    cout << "exact_sol_shift: "<<exact_sol_shift<<endl;
    cout << "T_period       : "<<T_period<<endl;
    printf("actual_end_time:%1.2f",simdata_->t_end_);
    cout <<"\nMax_iter: "<<simdata_->maxIter_<<endl;

    cout << "\nNumber of nodes: "<< grid_->Nfaces<<"  dx:  "<<grid_->dx<<endl;
    cout << "Scheme  order : "<< simdata_->scheme_order_  << endl;
    cout << "Runge-Kutta order : "<< simdata_->RK_order_    << endl;
    cout <<"===============================================\n";

    return;
}

void FDSolverAdvecDiffus::Reset_solver(){

    emptyarray(Nfaces_tot,Qn);
    emptyarray(grid_->N_exact_ppts,Q_exact_pp);
    emptyarray(Nfaces,Q_init);
    emptyarray(Nfaces,Q_exact);

    emptyarray(stencil_index);
    emptyarray(FD_coeff);

    emptyarray(stencil_index_2nd);
    emptyarray(FD_coeff_2nd);
    simdata_->Reset();

    return;
}

// Solver functions
//-------------------------------------------

void FDSolverAdvecDiffus::setup_stencil(){

    stencil_index = new int[scheme_order_+1];
    FD_coeff = new double [scheme_order_+1];

    stencil_index_2nd = new int[scheme_order_+1];
    FD_coeff_2nd = new double [scheme_order_+1];

    if(scheme_order_==1){      // first order upwind scheme

        _notImplemented("There is no 1st order stencil defined for advec-diffus");
        FatalError_exit("stencil setup");

    } else if (scheme_order_==2) { // 2nd order central scheme

        stencil_index[0] =  1;
        stencil_index[1] =  0;
        stencil_index[2] = -1;           //[j+1,j,j-1];

        stencil_index_2nd[0] = stencil_index[0];
        stencil_index_2nd[1] = stencil_index[1];
        stencil_index_2nd[2] = stencil_index[2];

        FD_coeff [0] =  0.5;
        FD_coeff [1] =  0.0;
        FD_coeff [2] = -0.5;

        FD_coeff_2nd [0] =  1.0;
        FD_coeff_2nd [1] = -2.0;
        FD_coeff_2nd [2] =  1.0;

    }else if (scheme_order_==3){

        _notImplemented("There is no 3rd order stencil defined for advec-diffus");
        FatalError_exit("stencil setup");

    }else if (scheme_order_==4){ // 4th order central scheme

        stencil_index[0] =  2;
        stencil_index[1] =  1;
        stencil_index[2] =  0;
        stencil_index[3] = -1;
        stencil_index[4] = -2;  //[j+2,j+1,j,j-1,j-2];

        stencil_index_2nd[0] = stencil_index[0];
        stencil_index_2nd[1] = stencil_index[1];
        stencil_index_2nd[2] = stencil_index[2];
        stencil_index_2nd[3] = stencil_index[3];
        stencil_index_2nd[4] = stencil_index[4];

        FD_coeff [0] =  -1.0/12.0;
        FD_coeff [1] =   2.0/3.0;
        FD_coeff [2] =   0.0;
        FD_coeff [3] =  -2.0/3.0;
        FD_coeff [4] =   1.0/12.0;

        FD_coeff_2nd [0] =  -1.0/12.0;
        FD_coeff_2nd [1] =   4.0/3.0;
        FD_coeff_2nd [2] =  -2.50;
        FD_coeff_2nd [3] =   4.0/3.0;
        FD_coeff_2nd [4] =  -1.0/12.0;

    }else if (scheme_order_==6){ // 6th order central scheme

        stencil_index[0] =  3;
        stencil_index[1] =  2;
        stencil_index[2] =  1;
        stencil_index[3] =  0;
        stencil_index[4] = -1;
        stencil_index[5] = -2;
        stencil_index[6] = -3;  //[j+3,j+2,j+1,j,j-1,j-2,j-3];

        stencil_index_2nd[0] = stencil_index[0];
        stencil_index_2nd[1] = stencil_index[1];
        stencil_index_2nd[2] = stencil_index[2];
        stencil_index_2nd[3] = stencil_index[3];
        stencil_index_2nd[4] = stencil_index[4];
        stencil_index_2nd[5] = stencil_index[5];
        stencil_index_2nd[6] = stencil_index[6];

        FD_coeff [0] =   1.0/6.0;
        FD_coeff [1] =  -0.15;
        FD_coeff [2] =   0.75;
        FD_coeff [3] =   0.00;
        FD_coeff [4] =  -0.75;
        FD_coeff [5] =   0.15;
        FD_coeff [6] =  -1.0/6.0;

        FD_coeff_2nd [0] =   1.0/90.0;
        FD_coeff_2nd [1] =   0.15;
        FD_coeff_2nd [2] =   1.50;
        FD_coeff_2nd [3] =  -49.0/18.0;
        FD_coeff_2nd [4] =   1.50;
        FD_coeff_2nd [5] =   0.15;
        FD_coeff_2nd [6] =   1.0/90.0;
    }

    return;
}

void FDSolverAdvecDiffus::CalcTimeStep(){

    T_period = (grid_->xf - grid_->x0) / simdata_->a_wave_;

    if(simdata_->calc_dt_flag==1){

        CFL = simdata_->CFL_;
        time_step = (grid_->dx * CFL )/ simdata_->a_wave_;
        last_time_step = time_step;
        simdata_->dt_ = time_step;

    }else if(simdata_->calc_dt_flag==0){

        time_step = simdata_->dt_;
        last_time_step = time_step;
        CFL = simdata_->a_wave_ * time_step / grid_->dx ;
        simdata_->CFL_ = CFL;

    }else {

        FatalError_exit("Wrong Calc_dt_flag");
    }

    // Determining end of simulation parameters:
    //----------------------------------------------------

    if(simdata_->end_of_sim_flag_==0){

        simdata_->t_end_ = simdata_->Nperiods * T_period;

        simdata_->maxIter_ = (int) ceil(simdata_->t_end_/time_step);

        if((simdata_->maxIter_ * time_step) > simdata_->t_end_ ){

            last_time_step = simdata_->t_end_ - ((simdata_->maxIter_-1) * time_step);

        }else if((simdata_->maxIter_ * time_step) < (simdata_->Nperiods * T_period) ){

            last_time_step = simdata_->t_end_ - (simdata_->maxIter_ * time_step);
        }

    }else if(simdata_->end_of_sim_flag_==1){

        simdata_->Nperiods = simdata_->t_end_/T_period;
        simdata_->maxIter_ = (int) ceil(simdata_->t_end_/time_step);

        if((simdata_->maxIter_ * time_step) > simdata_->t_end_ ){

            last_time_step = simdata_->t_end_ - ((simdata_->maxIter_-1) * time_step);

        }else if((simdata_->maxIter_ * time_step) < (simdata_->Nperiods * T_period) ){

            last_time_step = simdata_->t_end_ - (simdata_->maxIter_ * time_step);
        }

    }else if(simdata_->end_of_sim_flag_==2){

        simdata_->t_end_ = simdata_->maxIter_ * time_step;
        simdata_->Nperiods = simdata_->t_end_/T_period;

    }else{
        FatalError_exit("Wrong end_of_simulation_flag");
    }

    return;
}

void FDSolverAdvecDiffus::InitSol(){

    register int j;

    int k=0;

    for(j=0; j<Nfaces; j++){

        for(k=0; k<Ndof; k++){

            if(simdata_->wave_form_==3)
                Q_init[j][k] =eval_init_u_decay_burger_turb(grid_->X[j]);
            else
                Q_init[j][k] = eval_init_sol(grid_->X[j]);

            Qn[j+Nghost_l][k] = Q_init[j][k];
        }
    }

    return;
}

void FDSolverAdvecDiffus::ComputeExactSolShift(){

    // Preparing shift information:
    //-------------------------------
    double a=0.;

    wave_length_ = grid_->xf - grid_->x0 ;
    a = simdata_->a_wave_;
    exact_sol_shift = (a * simdata_->t_end_ );

    return;
}

void FDSolverAdvecDiffus::update_ghost_sol(double **Qn_){

    register int i;

    int k;

    // updating the left boundary:
    for(i=0; i<Nghost_l; i++)
        for(k=0; k<Ndof; k++)
            Qn_[i][k] = Qn_[Nfaces-1+i][k]; //(Nghost_l+Nfaces-1)-Nghost_l+i

    //updating the right boundary:
    for(i=0; i<Nghost_r; i++)
        for(k=0; k<Ndof; k++)
            Qn_[i+Nghost_l+Nfaces][k] = Qn_[i+Nghost_l+1][k];

    return;
}

void FDSolverAdvecDiffus::UpdateResid(double **Resid_, double **Qn_){

    // First Update ghost nodes:
    //-----------------------------

    update_ghost_sol(Qn_);

    // Nodes loop to calculate and update the residual:
    //----------------------------------------------------

    register int i;

    int k=0,j=0,s1,s2;

    double Idx = grid_->Idx;

    double temp_inv=0.0, temp_visc=0.0, invFlux=0.0;

    for(i=0; i<Nfaces; i++){
        for(k=0; k<Ndof; k++){
            temp_inv=0.0;
            temp_visc=0.0;
            for(j=0; j<scheme_order_+1; j++){
                s1 = stencil_index[j];
                s2 = stencil_index_2nd[j];
                invFlux = evaluate_inviscid_flux(Qn_[i+Nghost_l+s1][k]);
                temp_inv  += invFlux * FD_coeff[j];
                temp_visc += Qn_[i+Nghost_l+s2][k] * FD_coeff_2nd[j];
            }

            Resid_[i][k] = (- temp_inv * Idx )
                    + ( temp_visc *Idx*Idx *simdata_->thermal_diffus );
        }
    }

    return;
}

double FDSolverAdvecDiffus::evaluate_inviscid_flux(const double& qn_){

    if(simdata_->wave_form_==2 || simdata_->wave_form_==3){ // Burgers equation
        return (0.5 *qn_*qn_);
    }else if(simdata_->wave_form_==0 || simdata_->wave_form_==1){
        return (qn_*simdata_->a_wave_);
    }else{
        return 0.0;
    }
}

void FDSolverAdvecDiffus::Compute_exact_sol_for_plot(){

    register int j;

    double xx=0.0;
    double x0,x1;

    if(simdata_->wave_form_==1){

        for(j=0; j<grid_->N_exact_ppts; j++){

            xx = grid_->x_exact_ppts[j]- exact_sol_shift;

            x0 = xx - wave_length_*floor(xx/wave_length_);
            x1 = xx + wave_length_*floor(xx/-wave_length_);

            if(x0==0 && x1==0)
                Q_exact_pp[j][0] = 0.5*(eval_init_sol(x0)+ eval_init_sol(x1));
            else
                Q_exact_pp[j][0] = (eval_init_sol(x0)+ eval_init_sol(x1));
        }

    }else if(simdata_->wave_form_==0){

       for(j=0; j<grid_->N_exact_ppts; j++){

            xx = grid_->x_exact_ppts[j]- exact_sol_shift;
            Q_exact_pp[j][0] = eval_init_sol(xx);
       }

    }else if(simdata_->wave_form_==3){

       for(j=0; j<grid_->N_exact_ppts; j++){

            xx = grid_->x_exact_ppts[j]- exact_sol_shift;
            Q_exact_pp[j][0] = eval_init_u_decay_burger_turb(xx);
       }

    }

    return;
}

void FDSolverAdvecDiffus::Compute_exact_sol(){

    register int j;

    double xx=0.0;
    double x0,x1;

    if(simdata_->wave_form_==1){

        for(j=0; j<Nfaces; j++){

            xx = grid_->X[j] - exact_sol_shift;

            x0 = xx - wave_length_*floor(xx/wave_length_);
            x1 = xx + wave_length_*floor(xx/-wave_length_);

            if(x0==0 && x1==0)
                Q_exact[j][0] = 0.5*(eval_init_sol(x0)+ eval_init_sol(x1));
            else
                Q_exact[j][0] = (eval_init_sol(x0)+ eval_init_sol(x1));
        }

    }else if(simdata_->wave_form_==0){

       for(j=0; j<Nfaces; j++){

            xx = grid_->X[j]- exact_sol_shift;
            Q_exact[j][0] = eval_init_sol(xx);
       }

    }else if(simdata_->wave_form_==3){

       for(j=0; j<Nfaces; j++){

            xx = grid_->X[j]- exact_sol_shift;
            Q_exact[j][0] = eval_init_u_decay_burger_turb(xx);
       }

    }

    return;

}

double FDSolverAdvecDiffus::eval_init_sol(const double& xx){

    if(simdata_->wave_form_==0){

        return sin(simdata_->wave_freq_*PI*xx);  // wave_freq_ x PI

    }else if(simdata_->wave_form_==1){

        return exp(-simdata_->Gaussian_exponent_*pow(xx,2));
    }else{
        _notImplemented("Wave form is not implemented");
    }
}

double FDSolverAdvecDiffus::eval_init_u_decay_burger_turb(const double& xx_){

    register int i;
    double u_=0.;
    double dk_,k_max_, E_, k_, epsi_;
    k_max_ = simdata_->max_wave_no_;
    dk_ = 1.0;
    int n_pts_=k_max_/dk_;

    for(i=0; i<n_pts_; i++){
        k_ = simdata_->k_wave_no_[i];
        epsi_= simdata_->epsi_phase_[i];
        E_= simdata_->energy_spect_[i];
        u_ += sqrt(2.*E_ ) * cos (k_ * xx_ + 2.*PI*epsi_) ;
    }

    return u_;
}

double FDSolverAdvecDiffus::L1_error_nodal_sol(){

    register int j;

    double L1_error=0.0;

    for(j=0; j<Nfaces; j++)
        L1_error += fabs(Q_exact[j][0] - Qn[j+Nghost_l][0]);

    L1_error = L1_error/Nfaces;

    return L1_error;
}

double FDSolverAdvecDiffus::L2_error_nodal_sol(){

    register int j;

    double L2_error=0.0;

    for(j=0; j<Nfaces; j++)
        L2_error += pow((Q_exact[j][0] - Qn[j+Nghost_l][0]),2);

    L2_error = sqrt(L2_error/Nfaces);

    return L2_error;
}

void FDSolverAdvecDiffus::print_cont_vertex_sol(){

    register int j=0;

    char *fname=nullptr;
    fname = new char[150];

    if(simdata_->Sim_mode=="dt_const"){

        sprintf(fname,"%snodal/u_num_N%d_dt%1.3e_%1.3fT.dat"
                ,simdata_->case_postproc_dir
                ,grid_->Nelem
                ,time_step
                ,simdata_->Nperiods);

        FILE* sol_out=fopen(fname,"w");

        for(j=0; j<Nfaces; j++)
            fprintf(sol_out, "%2.10e %2.10e\n", grid_->X[j], Qn[j+Nghost_l][0]);

        fclose(sol_out);
        emptyarray(fname);

    }else{
        sprintf(fname,"%snodal/u_num_N%d_CFL%1.3f_%1.3fT.dat"
                ,simdata_->case_postproc_dir
                ,grid_->Nelem
                ,CFL
                ,simdata_->Nperiods);

        FILE* sol_out=fopen(fname,"w");

        for(j=0; j<Nfaces; j++)
            fprintf(sol_out, "%2.10e %2.10e\n"
                    ,grid_->X[j], Qn[j+Nghost_l][0]);

        fclose(sol_out);

        emptyarray(fname);
    }

    fname = new char[150];

    sprintf(fname,"%snodal/u_exact_%1.3fT.dat"
            ,simdata_->case_postproc_dir
            ,simdata_->Nperiods);

    FILE* sol_out=fopen(fname,"w");

    for(j=0; j<grid_->N_exact_ppts; j++)
        fprintf(sol_out, "%2.10e %2.10e\n"
                ,grid_->x_exact_ppts[j], Q_exact_pp[j][0]);

    fclose(sol_out);
    emptyarray(fname);

    return;
}

void FDSolverAdvecDiffus::dump_timeaccurate_sol(){

    register int j=0;

    char *fname=nullptr;
    fname = new char[400];
    sprintf(fname,"%stime_data/u_num_N%d_dt%1.3e_%1.3ft.dat"
            ,simdata_->case_postproc_dir
            ,grid_->Nelem
            ,time_step
            ,phy_time);

    FILE* sol_out=fopen(fname,"w");

    for(j=0; j<Nfaces; j++)
        fprintf(sol_out, "%2.10e %2.10e\n", grid_->X[j], Qn[j+Nghost_l][0]);

    fclose(sol_out);
    emptyarray(fname);

    return;
}

void FDSolverAdvecDiffus::dump_errors(double &L1_error_, double &L2_error_){

    char *fname=nullptr;
    fname = new char[150];

    if(simdata_->Sim_mode=="CFL_const"){

        sprintf(fname,"%serrors/nodal_errors_CFL%1.3f_%1.3fT.dat"
                ,simdata_->case_postproc_dir
                ,CFL
                ,simdata_->Nperiods);

        FILE* solerror_out=fopen(fname,"at+");

        fprintf(solerror_out, "%d %2.10e %2.10e\n"
                ,grid_->Nelem, L1_error_, L2_error_);

        fclose(solerror_out);
        emptyarray(fname); fname = new char[100];

        sprintf(fname,"%serrors/nodal_errors_N%d_CFL%1.3f_%1.3fT.dat"
                ,simdata_->case_postproc_dir
                ,grid_->Nelem
                ,CFL
                ,simdata_->Nperiods);

        solerror_out=fopen(fname,"w");

        fprintf(solerror_out,"%2.10e %2.10e\n",L1_error_, L2_error_);

        fclose(solerror_out);

    }else if(simdata_->Sim_mode=="dt_const"){

        sprintf(fname,"%serrors/nodal_errors_dt%1.3e_%1.3fT.dat"
                ,simdata_->case_postproc_dir
                ,time_step
                ,simdata_->Nperiods);

        FILE* solerror_out=fopen(fname,"at+");

        fprintf(solerror_out, "%d %2.10e %2.10e\n"
                ,grid_->Nelem, L1_error_, L2_error_);

        fclose(solerror_out);
        emptyarray(fname); fname = new char[100];

        sprintf(fname,"%serrors/nodal_errors_N%d_dt%1.3e_%1.3fT.dat"
                ,simdata_->case_postproc_dir
                ,grid_->Nelem
                ,time_step
                ,simdata_->Nperiods);

        solerror_out=fopen(fname,"w");

        fprintf(solerror_out,"%2.10e %2.10e\n",L1_error_, L2_error_);

        fclose(solerror_out);

    }else{

        sprintf(fname,"%serrors/nodal_errors_N%d_CFL%1.3f_%1.3fT.dat"
                ,simdata_->case_postproc_dir
                ,grid_->Nelem
                ,CFL
                ,simdata_->Nperiods);

        FILE* solerror_out=fopen(fname,"w");

        fprintf(solerror_out,"%2.10e %2.10e\n",L1_error_, L2_error_);

        fclose(solerror_out);
    }

    emptyarray(fname);

    return;
}


























