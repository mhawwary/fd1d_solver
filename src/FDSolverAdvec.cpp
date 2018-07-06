#include "FDSolverAdvec.hpp"

// Constructor/Destructor/ Setup functions:
//------------------------------------------------
FDSolverAdvec::FDSolverAdvec(void){}

FDSolverAdvec::~FDSolverAdvec(void){

    Reset_solver();
}

void FDSolverAdvec::setup_solver(GridData* meshdata_, SimData& osimdata_){

    grid_ = meshdata_;
    simdata_ = &osimdata_;

    Ndof= 1;

    Nfaces = grid_->Nfaces;
    n_linsys = Nfaces-1;

    scheme_type_ = simdata_->scheme_type_;
    scheme_order_ = simdata_->scheme_order_;
    filter_type_ = simdata_->filter_type_;
    filter_order_ = simdata_->filter_order_;
    filter_stencil_size_ = simdata_->filter_stencil_size_;
    filter_activate_flag = simdata_->filter_activate_flag_;
    filter_alpha_ = simdata_->filter_alpha_;

    if(scheme_type_=="explicit"){
        if(scheme_order_==1){
            Nghost_l = 1;
            Nghost_r=0;
            stencil_width_=2;

        }else if(scheme_order_==2){
            if(simdata_->upwind_biased_==0){
            Nghost_l=1;
            Nghost_r=1;
            }else if(simdata_->upwind_biased_==2){
                Nghost_l=2;
                Nghost_r=1;
            }
            stencil_width_=3;

        }else if(scheme_order_==3){
            stencil_width_=4;
            if(simdata_->upwind_biased_==1){ // 1-point biased
                Nghost_l=2;
                Nghost_r=1;
            }else if(simdata_->upwind_biased_==0){ // fully-upwind
                Nghost_l=3;
                Nghost_r=0;
            }

        }else if(scheme_order_==4){
            stencil_width_=5;
            if(simdata_->upwind_biased_==2){ // 2-point biased
                Nghost_l=3;
                Nghost_r=1;
            }else if(simdata_->upwind_biased_==0){ // central
                Nghost_l=2;
                Nghost_r=2;
            }

        }else if(scheme_order_==6){
            stencil_width_=7;
            if(simdata_->upwind_biased_==2){ // 2-point biased
                Nghost_l=4;
                Nghost_r=2;
            }else if(simdata_->upwind_biased_==0){ // central
                Nghost_l=3;
                Nghost_r=3;
            }
        }

    }else if(scheme_type_=="DRP4s7"
             || scheme_type_=="Rem2s7"){ //Tamwebb1993 & LindersNordstrom2015, central schemes
        Nghost_l=3;
        Nghost_r=3;
        stencil_width_=7;
        if(scheme_type_=="DRP4s7")
            scheme_order_=4;
        else
            scheme_order_=2;
    }else if(scheme_type_=="BB4s9"
             || scheme_type_=="Rem2s9"){ //BogeyBailly2004 & LindersNordstrom2015, central schemes
        Nghost_l=4;
        Nghost_r=4;
        stencil_width_=9;
        if(scheme_type_=="BB4s9")
            scheme_order_=4;
        else
            scheme_order_=2;
    }else if(scheme_type_=="BB4s11"){ //BogeyBailly2004, central schemes
        Nghost_l=5;
        Nghost_r=5;
        stencil_width_=11;
        scheme_order_=4;
    }else if(scheme_type_=="BB4s13"){ //BogeyBailly2004, central schemes
        Nghost_l=6;
        Nghost_r=6;
        stencil_width_=13;
        scheme_order_=4;

    }else if(scheme_type_=="implicit"){
        if(scheme_order_==4){
            Nghost_l=1;
            Nghost_r=1;

        }else if(scheme_order_==6){
            Nghost_l=2;
            Nghost_r=2;
        }
    }

    Nfaces_tot = Nghost_l + Nfaces + Nghost_r;

    if(filter_activate_flag==1){
        std::string Bound_type = "Periodic";
        if(simdata_->filter_type_=="pade"){
            filter = new PadeFilter();
            filter->setup_filter(simdata_->filter_type_,filter_order_,Nfaces
                                 ,filter_alpha_,Bound_type);
        }else if(simdata_->filter_type_=="BogeyBailly"
                 || simdata_->filter_type_=="standard"){
            filter = new ExplicitFilter();
            filter->setup_filter(simdata_->filter_type_,filter_stencil_size_,Nfaces
                                 ,filter_alpha_,Bound_type);
        }
    }

    Qn      =  new double* [Nfaces_tot];
    Q_init  =  new double* [Nfaces];
    Q_exact =  new double* [Nfaces];
    Q_exact_pp = new double*[grid_->N_exact_ppts];
    //qq_exact_time = new double*[Nfaces];

    register int i;

    for(i=0; i<Nfaces_tot; i++)
        Qn[i]     = new double[Ndof];

    for(i=0; i<Nfaces; i++){
        Q_init[i] = new double[Ndof];
        Q_exact[i] = new double[Ndof];
    }

    for(i=0; i<grid_->N_exact_ppts; i++)
        Q_exact_pp[i] = new double [Ndof];

    setup_coefficients();

    SetPhyTime(simdata_->t_init_);

    //Wave data:
    wave_length_ = simdata_->xf_-simdata_->x0_;
    wave_speed_ = simdata_->a_wave_;
    ComputeExactSolShift();

    //Compute exact solutions at time 0
    Compute_exact_sol();
    Compute_exact_sol_for_plot();

    return;
}

void FDSolverAdvec::Reset_solver(){

    emptyarray(Nfaces_tot,Qn);
    emptyarray(grid_->N_exact_ppts,Q_exact_pp);
    emptyarray(Nfaces,Q_init);
    emptyarray(Nfaces,Q_exact);
    //emptyarray(Nfaces,qq_exact_time);

    emptyarray(stencil_index);
    emptyarray(FD_coeff);

    simdata_->Reset();

    emptyarray(dfdx_);
    emptyarray(alpha_vec_f1_);
    emptyarray(b_vec_);
    emptyarray(RHS_f1_);

    emptypointer(filter);

    return;
}

// Solver functions
//-------------------------------------------

void FDSolverAdvec::setup_coefficients(){

    if(scheme_type_ == "explicit"){
        stencil_index = new int[stencil_width_];
        FD_coeff = new double [stencil_width_];

        if(scheme_order_==1){      // first order upwind scheme
            stencil_index[0] =  0;
            stencil_index[1] = -1;           //[j,j-1];
            FD_coeff [0] =  1.0;
            FD_coeff [1] = -1.0;

        } else if (scheme_order_==2) { // 2nd order central scheme
            if(simdata_->upwind_biased_==0) {  // central:
                stencil_index[0] =  1;
                stencil_index[1] =  0;
                stencil_index[2] = -1;           //[j+1,j,j-1];
                FD_coeff [0] =  0.5;
                FD_coeff [1] =  0.0;
                FD_coeff [2] = -0.5;
            }else if(simdata_->upwind_biased_==2) {  // fully upwind:
                stencil_index[0] =  0;
                stencil_index[1] = -1;
                stencil_index[2] = -2;           //[j,j-1,j-2];
                FD_coeff [0] =  1.5;
                FD_coeff [1] = -2.0;
                FD_coeff [2] = -0.5;
            }

        }else if (scheme_order_==3){

            if(simdata_->upwind_biased_==0) {  // fully upwind:
                stencil_index[0] =  0;
                stencil_index[1] = -1;
                stencil_index[2] = -2;
                stencil_index[3] = -3;  // [j,j-1,j-2,j-3];
                FD_coeff [0] =  11.0/6.0;
                FD_coeff [1] =  -3.0;
                FD_coeff [2] =   3.0/2.0;
                FD_coeff [3] =  -1.0/3.0;

            }else if(simdata_->upwind_biased_==1){ // 1-point biased :
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

        }else if (scheme_order_==6){ // 6th order central scheme

            if(simdata_->upwind_biased_==0) {  // central
                stencil_index[0] =  3;
                stencil_index[1] =  2;
                stencil_index[2] =  1;
                stencil_index[3] =  0;
                stencil_index[4] = -1;
                stencil_index[5] = -2;
                stencil_index[6] = -3;  //[j+3,j+2,j+1,j,j-1,j-2,j-3];

                FD_coeff [0] =   1.0/60.0;
                FD_coeff [1] =  -0.15;
                FD_coeff [2] =   0.75;
                FD_coeff [3] =   0.00;
                FD_coeff [4] =  -0.75;
                FD_coeff [5] =   0.15;
                FD_coeff [6] =  -1.0/60.0;
            }else if(simdata_->upwind_biased_==2) {  // 2-point biased
                stencil_index[0] =  2;
                stencil_index[1] =  1;
                stencil_index[2] =  0;
                stencil_index[3] = -1;
                stencil_index[4] = -2;
                stencil_index[5] = -3;
                stencil_index[6] = -4;  //[j+2,j+1,j,j-1,j-2,j-3,j-4];

                FD_coeff [0] =  -1.0/30.0;
                FD_coeff [1] =  0.40;
                FD_coeff [2] =  7.0/12.0;
                FD_coeff [3] =  -4.0/3.0;
                FD_coeff [4] =  0.50;
                FD_coeff [5] =  -2.0/15.0;
                FD_coeff [6] =  1.0/60.0;
            }
        }

    }else if(scheme_type_=="DRP4s7"
             || scheme_type_=="Rem2s7"){
        stencil_index = new int[stencil_width_];
        FD_coeff = new double [stencil_width_];
        stencil_index[0] =  3;
        stencil_index[1] =  2;
        stencil_index[2] =  1;
        stencil_index[3] =  0;
        stencil_index[4] = -1;
        stencil_index[5] = -2;
        stencil_index[6] = -3;  //[j+3,j+2,j+1,j,j-1,j-2,j-3];

        if(scheme_order_==4||scheme_type_=="DRP4s7"){ //DRP4s7 TamWebb1993
            //Coefficient are from Cunha2014, Sjorgeen2017 optimized for K[0,1.1]
            FD_coeff [0] =   0.77088238051822552;
            FD_coeff [1] =  -0.166705904414580469;
            FD_coeff [2] =   0.02084314277031176;
            FD_coeff [3] =   0.00;
            FD_coeff [4] =  -0.02084314277031176;
            FD_coeff [5] =   0.166705904414580469;
            FD_coeff [6] =  -0.77088238051822552;
        }else if(scheme_order_==2||scheme_type_=="Rem2s7"){ //Rem2s7, LindersNordstrom2015
            FD_coeff [0] =   0.78028389;
            FD_coeff [1] =  -0.17585010;
            FD_coeff [2] =   0.02380544;
            FD_coeff [3] =   0.00;
            FD_coeff [4] =  -0.02380544;
            FD_coeff [5] =   0.17585010;
            FD_coeff [6] =  -0.78028389;
        }

    }else if(scheme_type_ == "implicit"){
        register int i;
        n_linsys = Nfaces-1;
        dfdx_   =  new double[n_linsys];
        alpha_vec_f1_ = new double[n_linsys];
        b_vec_ = new double[n_linsys];
        RHS_f1_ =  new double[n_linsys];

        stencil_index = new int[scheme_order_-2];

        if(scheme_order_==4){
            alpha_f1_ = 0.25;
            a_f1_ = 1.5 * 0.5;   // a * 0.5
            b_f1_ = 0.0;   // 0.0

            stencil_index[0] =  1;
            stencil_index[1] = -1;       //[j+1,j-1];

        }else if(scheme_order_==6){
            alpha_f1_ = 1.0/3.0;
            a_f1_ = 0.50 * 14.0/9.0;     // a * 0.5
            b_f1_ = 0.25 * 1.0/9.0;      // b * 0.25

            stencil_index[0] =  2;
            stencil_index[1] =  1;
            stencil_index[2] = -1;
            stencil_index[3] = -2;  //[j+2,j+1,j-1,j-2];
        }

        for(i=0; i<n_linsys; i++){
            alpha_vec_f1_[i] = alpha_f1_;
            b_vec_[i] = 1.0;
            RHS_f1_[i]=0.0;
        }

    }else{
        FatalError_exit("Wrong Scheme type for space solver\
                        , use either explicit or implicit");
    }

    return;
}

void FDSolverAdvec::CalcTimeStep(){

    T_period = (grid_->xf - grid_->x0) / simdata_->a_wave_;

    if(simdata_->calc_dt_flag==1){  // sepecify CFL

        CFL = simdata_->CFL_;
        time_step = (grid_->dx * CFL )/ simdata_->a_wave_;
        last_time_step = time_step;
        simdata_->dt_ = time_step;

    }else if(simdata_->calc_dt_flag==0){  // sepcify dt

        time_step = simdata_->dt_;
        last_time_step = time_step;
        CFL = simdata_->a_wave_ * time_step / grid_->dx ;
        simdata_->CFL_ = CFL;

    }else {

        FatalError_exit("Wrong Calc_dt_flag");
    }

    // Determining end of simulation parameters:
    //----------------------------------------------------
    double temp_tol=1e-8;
    if(simdata_->end_of_sim_flag_==0){ // use no. of periods

        simdata_->t_end_ = simdata_->Nperiods * T_period;

        simdata_->maxIter_ = (int) floor(simdata_->t_end_/time_step);

        if((simdata_->maxIter_ * time_step)
                > (simdata_->Nperiods * T_period) ){

            last_time_step = simdata_->t_end_
                    - ((simdata_->maxIter_-1) * time_step);

        }else if((simdata_->maxIter_ * time_step)
                 < (simdata_->Nperiods * T_period) ){

            last_time_step = simdata_->t_end_
                    - (simdata_->maxIter_ * time_step);
        }

    }else if(simdata_->end_of_sim_flag_==1){  // use final time

        simdata_->Nperiods = simdata_->t_end_/T_period;
        simdata_->maxIter_ = (int) floor(simdata_->t_end_/time_step);

        if((simdata_->maxIter_ * time_step)
                > (simdata_->t_end_-temp_tol) ){

            last_time_step = simdata_->t_end_ - ((simdata_->maxIter_-1) * time_step);

        }else if((simdata_->maxIter_ * time_step)
                 < (simdata_->t_end_+temp_tol) ){

            last_time_step = simdata_->t_end_ - (simdata_->maxIter_ * time_step);
        }

        if(last_time_step<=1e-10){
            last_time_step=time_step;
            simdata_->maxIter_--;
        }

    }else if(simdata_->end_of_sim_flag_==2){  // use no. of iterations

        simdata_->t_end_ = simdata_->maxIter_ * time_step;
        simdata_->Nperiods = simdata_->t_end_/T_period;

    }else{
        FatalError_exit("Wrong end_of_simulation_flag");
    }

    // Screen Output of input and simulation parameters:
    cout <<"\n===============================================\n";
    cout << "max eigenvalue : "<<max_eigen_advec<<endl;
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
    cout << "Scheme  order    : "<< simdata_->scheme_order_  << endl;
    cout << "Time scheme type : "<< simdata_->time_scheme_type_<<endl;
    if(simdata_->time_scheme_type_=="RungeKutta")
        cout << "Runge-Kutta order: "<< simdata_->RK_order_    << endl;
    cout <<"===============================================\n";

    return;
}

void FDSolverAdvec::InitSol(){

    register int j;
    int k=0;
    max_eigen_advec = simdata_->a_wave_;

    if(simdata_->eqn_type_=="linear_advec"){

        for(j=0; j<Nfaces; j++){
            for(k=0; k<Ndof; k++){
                Q_init[j][k] = eval_init_sol(grid_->X[j]);
                Qn[j+Nghost_l][k] = Q_init[j][k];
            }
        }
    }else if(simdata_->eqn_type_=="inv_burger"){
        max_eigen_advec = 0.0;
        for(j=0; j<Nfaces; j++){
            for(k=0; k<Ndof; k++){
                Q_init[j][k] = eval_init_sol(grid_->X[j]);
                Qn[j+Nghost_l][k] = Q_init[j][k];

                if(fabs(Q_init[j][k])>max_eigen_advec)
                    max_eigen_advec = fabs(Q_init[j][k]);
            }
        }
    }
    CalcTimeStep();

    return;
}

void FDSolverAdvec::ComputeExactSolShift(){
    exact_sol_shift = wave_speed_*phy_time;
    return;
}

void FDSolverAdvec::update_ghost_sol(double **Qn_){

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

void FDSolverAdvec::UpdateResid(double **Resid_, double **qn_){

    register int i;
    double Idx = grid_->Idx;  // 1/h
    int k=0;

    // First Update ghost nodes:
    //-----------------------------
    update_ghost_sol(qn_);
    // Nodes loop to calculate and update the residual:
    //----------------------------------------------------

    if(scheme_type_=="explicit"
            || scheme_type_=="DRP4s7"
            || scheme_type_=="Rem2s7"){
        // Nodes loop to calculate and update the residual:
        //----------------------------------------------------
        int j=0,s1;
        double temp_inv=0.0, invFlux=0.0;

        for(i=0; i<Nfaces; i++){
            for(k=0; k<Ndof; k++){
                temp_inv=0.0;
                for(j=0; j<stencil_width_; j++){
                    s1 = stencil_index[j];
                    invFlux = evaluate_inviscid_flux(qn_[i+Nghost_l+s1][k]);
                    temp_inv  += invFlux * FD_coeff[j];
                }
                Resid_[i][k] = - ( temp_inv * Idx );
            }
        }

    }else if(scheme_type_=="implicit"){
        // Nodes loop to calculate and update the residual:
        //----------------------------------------------------
        register int i;
        // compute 1st derivative:
        compute_RHS_f1_implicit(Idx, &qn_[Nghost_l], RHS_f1_);
        cyclic_tridiag_solve_mh(n_linsys,alpha_vec_f1_,b_vec_
                                ,alpha_vec_f1_,RHS_f1_,dfdx_);

        for(i=0; i<n_linsys; i++)  // n_linsys == Nfaces-1
            Resid_[i][0] = - dfdx_[i];

        Resid_[Nfaces-1][0] = Resid_[0][0];
    }

    return;
}

void FDSolverAdvec::filter_solution(double **qn_){
    filter->filtered_sol(&qn_[Nghost_l]);  // filtering the solution
    return;
}

void FDSolverAdvec::compute_RHS_f1_implicit(const double& Idx_, double** qn_
                                            , double*& RHS_temp_){
    // Calculation of the RHS for the f' equation:
    //-----------------------------------------------
    register int i;
    double fp1=0.0,fm1=0.0,fp2=0.0,fm2=0.0;

    // note that Idx_ = 1/dx

    for(i=0; i<n_linsys; i++){
        fp1 = evaluate_inviscid_flux(qn_[i+1][0]);
        fm1 = evaluate_inviscid_flux(qn_[i-1][0]);
        RHS_temp_[i] = a_f1_ *(fp1-fm1);
        if(scheme_order_==6){
            fp2 = evaluate_inviscid_flux(qn_[i+2][0]);
            fm2 = evaluate_inviscid_flux(qn_[i-2][0]);
            RHS_temp_[i] +=  b_f1_ *(fp2-fm2);
        }
        RHS_temp_[i] = RHS_temp_[i]*Idx_;
    }

    return;
}

double FDSolverAdvec::evaluate_inviscid_flux(const double& qn_){

    if(simdata_->eqn_type_=="inv_burger"){ // Burgers equation
        return (0.5 *qn_*qn_);
    }else if(simdata_->eqn_type_=="linear_advec"){ // linear_advection
        return (qn_*simdata_->a_wave_);
    }else{
        return 0.0;
    }
}

void FDSolverAdvec::Compute_exact_sol_for_plot(){

    register int j;

    double xx=0.0;
    double x0,x1;

    ComputeExactSolShift();

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

            xx = grid_->x_exact_ppts[j]-  exact_sol_shift;
            Q_exact_pp[j][0] = eval_init_sol(xx);
        }

    }

    return;
}

void FDSolverAdvec::Compute_exact_sol(){

    register int j;
    double xx=0.0;
    double x0,x1;

    ComputeExactSolShift();

    if(simdata_->wave_form_==1){  // Gaussian wave
        for(j=0; j<Nfaces; j++){
            xx = grid_->X[j] - exact_sol_shift;
            x0 = xx - wave_length_*floor(xx/wave_length_);
            x1 = xx + wave_length_*floor(xx/-wave_length_);

            if(x0==0 && x1==0)
                Q_exact[j][0] = 0.5*(eval_init_sol(x0)+ eval_init_sol(x1));
            else
                Q_exact[j][0] = (eval_init_sol(x0)+ eval_init_sol(x1));
        }
    }else if(simdata_->wave_form_==0){ // sine wave
        for(j=0; j<Nfaces; j++){
            xx = grid_->X[j]- exact_sol_shift;
            Q_exact[j][0] = eval_init_sol(xx);
        }
    }
    return;
}

void FDSolverAdvec::Compute_TimeAccurate_exact_sol(){

    register int j;
    double xx=0.0;
    double x0,x1;

    ComputeExactSolShift();

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

double FDSolverAdvec::eval_init_sol(const double& xx){

    if(simdata_->wave_form_==0){

        return sin(simdata_->wave_freq_*PI*xx);  // wave_freq_ x PI

    }else if(simdata_->wave_form_==1){

        return exp(-simdata_->Gaussian_exponent_*pow(xx,2));
    }else{
        _notImplemented("Wave form is not implemented");
    }
}

double FDSolverAdvec::L1_error_nodal_sol(){

    register int j;

    double L1_error=0.0;

    for(j=0; j<Nfaces; j++)
        L1_error += fabs(Q_exact[j][0] - Qn[j+Nghost_l][0]);

    L1_error = L1_error/Nfaces;

    return L1_error;
}

double FDSolverAdvec::L2_error_nodal_sol(){

    register int j;

    double L2_error=0.0;

    for(j=0; j<Nfaces; j++)
        L2_error += pow((Q_exact[j][0] - Qn[j+Nghost_l][0]),2);

    L2_error = sqrt(L2_error/Nfaces);

    return L2_error;
}

double FDSolverAdvec::dissipation_error(){

    register int j;

    double dissip_error=0.0, Qn_aver=0.0
            , Qex_aver=0.0, Qn_var=0.0, Qex_var=0.0;

    for(j=0; j<Nfaces; j++){
        Qn_aver  += Qn[j+Nghost_l][0];
        Qex_aver += Q_exact[j][0];
    }
    Qn_aver = Qn_aver/Nfaces;
    Qex_aver = Qex_aver/ Nfaces;

    for(j=0; j<Nfaces; j++){
        Qn_var += pow((Qn[j+Nghost_l][0]-Qn_aver),2);
        Qex_var += pow((Q_exact[j][0]-Qex_aver),2);
    }
    Qn_var = sqrt(Qn_var/Nfaces);
    Qex_var = sqrt(Qex_var/Nfaces);

    dissip_error = sqrt(pow((Qn_var-Qex_var),2) + pow((Qn_aver-Qex_aver),2));

    return dissip_error;
}

void FDSolverAdvec::print_cont_vertex_sol(){

    register int j=0;

    char *fname=nullptr;
    fname = new char[100];

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

    fname = new char[100];

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

void FDSolverAdvec::dump_timeaccurate_sol(){

    register int j=0;

    char *fname=nullptr;
    fname = new char[400];

    if(simdata_->Sim_mode=="normal"
            || simdata_->Sim_mode=="test"
            || simdata_->Sim_mode=="CFL_const"
            || simdata_->Sim_mode=="error_analysis_CFL"){
        sprintf(fname,"%stime_data/u_num_N%d_CFL%1.4f_%1.3ft.dat"
                ,simdata_->case_postproc_dir
                ,grid_->Nelem
                ,CFL
                ,phy_time);
    }else if(simdata_->Sim_mode=="dt_const"
             || simdata_->Sim_mode=="error_analysis_dt" ){
        sprintf(fname,"%stime_data/u_num_N%d_dt%1.3e_%1.3ft.dat"
                ,simdata_->case_postproc_dir
                ,grid_->Nelem
                ,time_step
                ,phy_time);
    }

    FILE* sol_out=fopen(fname,"w");

    for(j=0; j<Nfaces; j++)
        fprintf(sol_out, "%2.10e %2.10e\n", grid_->X[j], Qn[j+Nghost_l][0]);

    fclose(sol_out);
    emptyarray(fname);

    Compute_exact_sol_for_plot();
    fname = new char[400];
    sprintf(fname,"%stime_data/u_exact_%1.3ft.dat"
            ,simdata_->case_postproc_dir
            ,phy_time);

    FILE* sol_out1=fopen(fname,"w");

    for(j=0; j<grid_->N_exact_ppts; j++)
        fprintf(sol_out1, "%2.10e %2.10e\n"
                ,grid_->x_exact_ppts[j], Q_exact_pp[j][0]);

    fclose(sol_out1);
    emptyarray(fname);

    return;
}

void FDSolverAdvec::dump_timeaccurate_errors(){

    //Compute_TimeAccurate_exact_sol();
    Compute_exact_sol();
    double L1_error_=L1_error_nodal_sol();
    double L2_error_=L2_error_nodal_sol();

    char *fname=nullptr;
    fname = new char[100];
    sprintf(fname,"%serrors/errors_N%d_CFL%1.4f_%1.3fT.dat"
            ,simdata_->case_postproc_dir
            ,grid_->Nelem
            ,CFL
            ,simdata_->Nperiods);

    FILE* solerror_out=fopen(fname,"at+");

    fprintf(solerror_out, "%1.10f %2.10e %2.10e\n"
            ,phy_time, L1_error_, L2_error_);

    fclose(solerror_out);
    emptyarray(fname);

    printf("\tL1_nodal: %2.5e\t L2_nodal: %2.5e",L1_error_,L2_error_);

    return;
}

void FDSolverAdvec::dump_errors(double &L1_error_, double &L2_error_){

    char *fname=nullptr;
    fname = new char[100];

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






