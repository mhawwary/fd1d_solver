#include "ExplicitTimeSolver.hpp"

ExplicitTimeSolver::ExplicitTimeSolver(void){}

ExplicitTimeSolver::~ExplicitTimeSolver(void){

    Reset_time_solver();

}

void ExplicitTimeSolver::setupTimeSolver(FDSolver *fd_solver_, SimData *osimdata_){

    space_solver = fd_solver_;

    simdata = osimdata_;

    Nfaces = space_solver->GetNfaces();

    Nghost_l = space_solver->GetNghost_l();

    Ndof = space_solver->GetNdof();

    resid = new double* [Nfaces];

    time_scheme_type = simdata->time_scheme_type_;

    register int i;

    for(i=0; i<Nfaces; i++)
        resid[i] = new double[Ndof];

    dt_ = space_solver->GetTimeStep();

    if(time_scheme_type=="RungeKutta"){
        if(simdata->RK_order_>1){

            q_temp = new double*[Nfaces];

            for(i=0; i<Nfaces; i++)
                q_temp[i] = new double[Ndof];

            if(simdata->RK_order_==4){

                resid_temp0 = new double* [Nfaces];
                resid_temp1 = new double* [Nfaces];
                resid_temp2 = new double* [Nfaces];

                for(i=0; i<Nfaces; i++){
                    resid_temp0[i] = new double[Ndof];
                    resid_temp1[i] = new double[Ndof];
                    resid_temp2[i] = new double[Ndof];
                }
            }
        }

    }else if(time_scheme_type=="leapfrog"){
        q_temp = new double*[Nfaces];
        qnm1_temp = new double*[Nfaces];

        for(i=0; i<Nfaces; i++){
            q_temp[i] = new double[Ndof];
            qnm1_temp[i] = new double[Ndof];
        }
    }else if(time_scheme_type=="upwind_leapfrog"){
        q_temp = new double*[Nfaces];
        qnm1_temp = new double*[Nfaces+Nghost_l];

        for(i=0; i<Nfaces; i++)
            q_temp[i] = new double[Ndof];
        for(i=0; i<Nfaces+Nghost_l; i++)
            qnm1_temp[i] = new double[Ndof];
    }


    return;
}

void ExplicitTimeSolver::SolveOneStep(double **qn_){

    if(time_scheme_type=="RungeKutta"){
        switch (simdata->RK_order_) {

        case 1:

            FwdEuler(qn_);

            break;

        case 2:

            SSPRK22(qn_);

            break;

        case 3:

            SSPRK33(qn_);

            break;

        case 4:

            classicRK4(qn_);

            break;

        default:
            char *ss=nullptr; ss= new char[100];
            sprintf(ss,"RK order of %d ",simdata->RK_order_);
            _notImplemented(ss);
            emptyarray(ss);
            break;
        }
    }else if(time_scheme_type=="leapfrog"){

        if(IterNo==0){  // Solving the first time step using RK3
            CopyOldSol(qnm1_temp,qn_);
            SSPRK33(qn_);
        }else{
            leapfrog(qn_);
        }
    }else if(time_scheme_type=="upwind_leapfrog"){
        if(IterNo==0){  // Solving the first time step using RK3
            CopyOldSolwithGhost(qnm1_temp,qn_);
            SSPRK33(qn_);
            emptyarray(Nfaces,q_temp);
            q_temp = new double*[Nfaces+Nghost_l];
            for(register int i=0; i<Nfaces+Nghost_l; i++)
                q_temp[i] = new double[Ndof];
        }else{
            upwind_leapfrog(qn_);
        }
    }

    IterNo++;

    //    if(IterNo==simdata->maxIter_-1)
    //        dt_ = space_solver->GetLastTimeStep();

    return;
}

void ExplicitTimeSolver::CopyOldSol(double **q_t_, double **qn_){

    register int j;

    int k;

    for(j=0; j<Nfaces; j++)
        for(k=0; k<Ndof; k++)
            q_t_[j][k] = qn_[j+Nghost_l][k];

    return;
}

void ExplicitTimeSolver::CopyOldestSol(double **q_t_, double **qn_){

    register int j;

    int k;

    for(j=0; j<Nfaces; j++)
        for(k=0; k<Ndof; k++)
            q_t_[j][k] = qn_[j][k];

    return;
}

void ExplicitTimeSolver::CopyOldSolwithGhost(double **q_t_, double **qn_){

    register int j;

    int k;

    for(j=0; j<Nfaces+Nghost_l; j++)
        for(k=0; k<Ndof; k++)
            q_t_[j][k] = qn_[j][k];

    return;
}

void ExplicitTimeSolver::CopyOldResid(double **resid_t_, double **old_resid_){

    register int j;

    int k;

    for(j=0; j<Nfaces; j++)
        for(k=0; k<Ndof; k++)
            resid_t_[j][k] = old_resid_[j][k];

    return;
}

void ExplicitTimeSolver::leapfrog(double **q_){

    register int j;
    int k;

    CopyOldSol(q_temp,q_); // copying the oldest level

    for(j=0; j<Nfaces; j++)
        for(k=0; k<Ndof; k++)
            q_[j+Nghost_l][k] = qnm1_temp[j][k] + 2.0*dt_*resid[j][k];

    CopyOldestSol(qnm1_temp,q_temp); // from qn to q_n-1

    if(simdata->filter_activate_flag_==1)
        space_solver->filter_solution(q_);
    space_solver->UpdateResid(resid,q_);

    return;
}

void ExplicitTimeSolver::upwind_leapfrog(double **q_){

    register int j;
    int k;

    CopyOldSolwithGhost(q_temp,q_); // copying the oldest level

    for(j=0; j<Nfaces; j++)
        for(k=0; k<Ndof; k++)
            q_[j+Nghost_l][k] = qnm1_temp[j+Nghost_l-1][k]
                    + (q_temp[j+Nghost_l][k] - q_temp[j+Nghost_l-1][k])
                    + 2.0*dt_*resid[j][k];

    CopyOldSolwithGhost(qnm1_temp,q_temp); // from qn to q_n-1

    if(simdata->filter_activate_flag_==1)
        space_solver->filter_solution(q_);
    space_solver->UpdateResid(resid,q_);

    return;
}

void ExplicitTimeSolver::FwdEuler(double **q_){

    register int j;
    int k;

    for(j=0; j<Nfaces; j++)
        for(k=0; k<Ndof; k++)
            q_[j+Nghost_l][k] = q_[j+Nghost_l][k] + dt_ * resid[j][k];

    if(simdata->filter_activate_flag_==1)
        space_solver->filter_solution(q_);
    space_solver->UpdateResid(resid,q_);

    return;
}

void ExplicitTimeSolver::SSPRK22(double **q_){

    register int j;

    int k;

    CopyOldSol(q_temp,q_);  // Copying level n solution and saving it

    // Step1:
    //-----------
    for(j=0; j<Nfaces; j++)
        for(k=0; k<Ndof; k++){
            q_[j+Nghost_l][k] = q_temp[j][k] + dt_ * resid[j][k];
        }

    space_solver->UpdateResid(resid,q_);

    // Step2:
    //------------
    for(j=0; j<Nfaces; j++)
        for(k=0; k<Ndof; k++) {
            q_[j+Nghost_l][k] = 0.5 * ( q_temp[j][k] +  q_[j+Nghost_l][k]
                    + dt_ * resid[j][k] );
        }

    if(simdata->filter_activate_flag_==1)
        space_solver->filter_solution(q_);
    space_solver->UpdateResid(resid,q_);

    return;
}

void ExplicitTimeSolver::SSPRK33(double **q_){

    register int j;

    int k;

    CopyOldSol(q_temp,q_);  // Copying level n solution and saving it


    // Step1:
    //-----------
    for(j=0; j<Nfaces; j++)
        for(k=0; k<Ndof; k++)
            q_[j+Nghost_l][k] = q_temp[j][k] + dt_ * resid[j][k];

    space_solver->UpdateResid(resid,q_);

    // Step2:
    //------------
    for(j=0; j<Nfaces; j++)
        for(k=0; k<Ndof; k++)
            q_[j+Nghost_l][k] =  (0.75 * q_temp[j][k] )
                    + 0.25 * ( q_[j+Nghost_l][k] + dt_ * resid[j][k] );

    space_solver->UpdateResid(resid,q_);

    // Step3:
    //--------------
    for(j=0; j<Nfaces; j++)
        for(k=0; k<Ndof; k++)
            q_[j+Nghost_l][k] =  ( q_temp[j][k]/3. )
                    + 2. * ( q_[j+Nghost_l][k] + dt_ * resid[j][k] ) / 3.;

    if(simdata->filter_activate_flag_==1)
        space_solver->filter_solution(q_);
    space_solver->UpdateResid(resid,q_);

    return;
}

void ExplicitTimeSolver::classicRK4(double **q_){

    register int j;

    int k;

    CopyOldSol(q_temp,q_);  // Copying level 0 solution and saving it
    CopyOldResid(resid_temp0,resid);  // Copying level 0 residual and saving it

    // Step1:
    //-----------
    for(j=0; j<Nfaces; j++)
        for(k=0; k<Ndof; k++)
            q_[j+Nghost_l][k] = q_temp[j][k] + 0.5 * dt_ * resid[j][k];

    space_solver->UpdateResid(resid,q_);

    CopyOldResid(resid_temp1,resid);  // Copying level 1 residual and saving it

    // Step2:
    //------------
    for(j=0; j<Nfaces; j++)
        for(k=0; k<Ndof; k++)
            q_[j+Nghost_l][k] = q_temp[j][k] + 0.5 * dt_ * resid[j][k];

    space_solver->UpdateResid(resid,q_);

    CopyOldResid(resid_temp2,resid);  // Copying level 2 residual and saving it

    // Step3:
    //--------------
    for(j=0; j<Nfaces; j++)
        for(k=0; k<Ndof; k++)
            q_[j+Nghost_l][k] = q_temp[j][k] + dt_ * resid[j][k];

    space_solver->UpdateResid(resid,q_);

    // Step4:
    //--------------
    for(j=0; j<Nfaces; j++)
        for(k=0; k<Ndof; k++)
            q_[j+Nghost_l][k] = q_temp[j][k]
                    + (dt_/6.0) * ( resid_temp0[j][k] + 2.0*resid_temp1[j][k]
                                    + 2.0*resid_temp2[j][k] + resid[j][k] );

    if(simdata->filter_activate_flag_==1)
        space_solver->filter_solution(q_);
    space_solver->UpdateResid(resid,q_);

    return;
}

void ExplicitTimeSolver::ComputeInitialResid(double **qn_){

    //space_solver->filter_solution(qn_);
    space_solver->UpdateResid(resid,qn_);

    return;
}

void ExplicitTimeSolver::Reset_time_solver(){

    emptyarray(Nfaces,resid);

    if(time_scheme_type=="RungeKutta"){
        emptyarray(Nfaces,q_temp);
    }else if(time_scheme_type=="leapfrog"){
        emptyarray(Nfaces,q_temp);
        emptyarray(Nfaces,qnm1_temp);
    }else if(time_scheme_type=="upwind_leapfrog"){
        emptyarray(Nfaces+Nghost_l,q_temp);
        emptyarray(Nfaces+Nghost_l,qnm1_temp);
    }

    emptyarray(Nfaces,resid_temp0);
    emptyarray(Nfaces,resid_temp1);
    emptyarray(Nfaces,resid_temp2);

    return;
}

