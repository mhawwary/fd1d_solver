#include "FDSolver.hpp"

// Constructor/Destructor/ Setup functions:
//------------------------------------------------
FDSolver::FDSolver(void){}

FDSolver::~FDSolver(void){

    Reset_solver();
}

void FDSolver::setup_solver(GridData* meshdata_, SimData& osimdata_){

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

void FDSolver::Reset_solver(){

    emptyarray(Nfaces_tot,Qn);
    emptyarray(grid_->N_exact_ppts,Q_exact_pp);
    emptyarray(Nfaces,Q_init);
    emptyarray(Nfaces,Q_exact);

    emptyarray(stencil_index);
    emptyarray(FD_coeff);

    return;
}

// Solver functions
//-------------------------------------------

void FDSolver::setup_stencil(){

    stencil_index = new int[scheme_order_+1];
    FD_coeff = new double [scheme_order_+1];

    if(scheme_order_==1){      // first order upwind scheme

        stencil_index[0] =  0;
        stencil_index[1] = -1;           //[j,j-1];

        FD_coeff [0] =  1;
        FD_coeff [1] = -1;

    } else if (scheme_order_==2) { // 2nd order central scheme

        stencil_index[0] =  1;
        stencil_index[1] =  0;
        stencil_index[2] = -1;           //[j+1,j,j-1];

        FD_coeff [0] =  0.5;
        FD_coeff [1] =  0.0;
        FD_coeff [2] = -0.5;

    }else if (scheme_order_==3){

        if(simdata_->upwind_biased_==0) {  //for 3rd order fully upwind:
            stencil_index[0] =  0;
            stencil_index[1] = -1;
            stencil_index[2] = -2;
            stencil_index[3] = -3;  // [j,j-1,j-2,j-3];

            FD_coeff [0] =  11.0/6.0;
            FD_coeff [1] =  -3.0;
            FD_coeff [2] =   3.0/2.0;
            FD_coeff [3] =  -1.0/3.0;

        }else if(simdata_->upwind_biased_==1){ // for 3rd order upwind biased :
            stencil_index[0] =  1;
            stencil_index[1] =  0;
            stencil_index[2] = -1;
            stencil_index[3] = -2;  //[j+1,j,j-1,j-2];

            FD_coeff [0] =  1.0/3.0;
            FD_coeff [1] =  0.5;
            FD_coeff [2] = -1.0;
            FD_coeff [3] =  1.0/6.0;

        }else{
            FatalError("Wrong biased upwind parameter for 3rd order scheme");
        }

    }else if (scheme_order_==4){ // 4th order central scheme

        stencil_index[0] =  2;
        stencil_index[1] =  1;
        stencil_index[2] =  0;
        stencil_index[3] = -1;
        stencil_index[4] = -2;  //[j+2,j+1,j,j-1,j-2];

        FD_coeff [0] =  -1.0/12.0;
        FD_coeff [1] =   2.0/3.0;
        FD_coeff [2] =   0.0;
        FD_coeff [3] =  -2.0/3.0;
        FD_coeff [4] =   1.0/12.0;
    }

    return;
}

void FDSolver::CalcTimeStep(){

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

void FDSolver::InitSol(){

    register int j;

    int k=0;

    for(j=0; j<Nfaces; j++){

        for(k=0; k<Ndof; k++){

            Q_init[j][k] = eval_init_sol(grid_->X[j]);

            Qn[j+Nghost_l][k] = Q_init[j][k];
        }
    }

    return;
}

void FDSolver::ComputeExactSolShift(){

    // Preparing shift information:
    //-------------------------------
    double a=0.;

    wave_length_ = grid_->xf - grid_->x0 ;
    a = simdata_->a_wave_;
    exact_sol_shift = (a * simdata_->t_end_ );

    return;
}

void FDSolver::update_ghost_sol(double **Qn_){

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

void FDSolver::UpdateResid(double **Resid_, double **Qn_){

    // First Update ghost nodes:
    //-----------------------------

    update_ghost_sol(Qn_);

    // Nodes loop to calculate and update the residual:
    //----------------------------------------------------

    register int i;

    int k=0,j=0,s;

    double Idx = grid_->Idx;

    double temp=0.0;

    for(i=0; i<Nfaces; i++){
        for(k=0; k<Ndof; k++){
            temp=0.0;
            for(j=0; j<scheme_order_+1; j++){
                s = stencil_index[j];
                temp += Qn_[i+Nghost_l+s][k] * FD_coeff[j];
            }

            Resid_[i][k] = - temp * Idx;
        }
    }

    return;
}

void FDSolver::Compute_exact_sol_for_plot(){

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

    }

    return;
}

void FDSolver::Compute_exact_sol(){

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

    }

    return;

}

double FDSolver::eval_init_sol(const double& xx){

    if(simdata_->wave_form_==0){

        return sin(2*PI*xx);

    }else if(simdata_->wave_form_==1){

        return exp(-simdata_->Gaussian_exponent_*pow(xx,2));
    }else{
        _notImplemented("Wave form is not implemented");
    }
}

double FDSolver::ComputeSolNodalError(){

    register int j;

    double L2_error=0.0;

    for(j=0; j<Nfaces; j++)
        L2_error += pow((Q_exact[j][0] - Qn[j+Nghost_l][0]),2);

    L2_error = sqrt(L2_error/Nfaces);

    return L2_error;
}

void FDSolver::print_cont_vertex_sol(){

    register int j=0;

    char *fname=nullptr;
    fname = new char[200];

    sprintf(fname,"%snodal/u_num_N%d_CFL%1.3f_%1.1fT.dat"
            ,simdata_->case_postproc_dir, simdata_->Nelem_
            ,simdata_->CFL_
            ,simdata_->Nperiods);

    FILE* sol_out=fopen(fname,"w");

    for(j=0; j<Nfaces; j++)
        fprintf(sol_out, "%2.10e %2.10e\n", grid_->X[j], Qn[j+Nghost_l][0]);

    fclose(sol_out);

    emptyarray(fname);

    fname = new char[200];

    sprintf(fname,"%snodal/u_exact_N%d_CFL%1.3f_%1.1fT.dat"
            ,simdata_->case_postproc_dir,simdata_->Nelem_
            ,simdata_->CFL_
            ,simdata_->Nperiods);

    FILE* sol_out1=fopen(fname,"w");

    for(j=0; j<grid_->N_exact_ppts; j++)
        fprintf(sol_out1, "%2.10e %2.10e\n"
                ,grid_->x_exact_ppts[j], Q_exact_pp[j][0]);

    fclose(sol_out1);

    emptyarray(fname);

    return;
}

void FDSolver::dump_errors(double &L2_error_){

    char *fname=nullptr;
    fname = new char[200];

    sprintf(fname,"%serrors/sol_errors_CFL%1.3f_%1.1fT.dat"
            ,simdata_->case_postproc_dir
            ,simdata_->CFL_
            ,simdata_->Nperiods);

    FILE* solerror_out=fopen(fname,"w");

    fprintf(solerror_out, "%d %2.10e\n"
            ,grid_->Nelem, L2_error_);

     fclose(solerror_out);

     emptyarray(fname);

     // Dumping all errors in one file as a function of beta:
     //--------------------------------------------------------
//     fname = new char[200];

//     sprintf(fname,"%serrors/sol_errors_N%d_CFL%1.2f_allBeta_%1.1fT.dat"
//             ,simdata_->case_postproc_dir
//             ,grid_->Nelem
//             ,simdata_->CFL_
//             ,simdata_->Nperiods);

//     solerror_out=fopen(fname,"at+");

//     fprintf(solerror_out, "%1.2f %2.10e %2.10e\n"
//             ,simdata_->upwind_param_, proj_sol_L2, aver_L2);

//      fclose(solerror_out);

//      emptyarray(fname);

    return;
}


























