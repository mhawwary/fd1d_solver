#include"ExplicitTimeSolver.hpp"

ExplicitTimeSolver::ExplicitTimeSolver(){

 return;
}

ExplicitTimeSolver::~ExplicitTimeSolver(){

    Reset_time_solver();

 return;
}

void ExplicitTimeSolver::setupTimeSolver(FDSolver *fd_solver_, SimData *osimdata_){

    space_solver = fd_solver_;

    simdata = osimdata_;

    Nfaces = space_solver->GetNfaces();

    Nghost_l = space_solver->GetNghost_l();

    Ndof = space_solver->GetNdof();

    resid = new double* [Nfaces];

    register int i;

    for(i=0; i<Nfaces; i++)
        resid[i] = new double[Ndof];

    dt_ = space_solver->GetTimeStep();

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

    return;
}

void ExplicitTimeSolver::SolveOneStep(double **qn_){

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

    IterNo++;

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

void ExplicitTimeSolver::CopyOldResid(double **resid_t_, double **old_resid_){

    register int j;

    int k;

    for(j=0; j<Nfaces; j++)
        for(k=0; k<Ndof; k++)
            resid_t_[j][k] = old_resid_[j][k];

    return;
}

void ExplicitTimeSolver::FwdEuler(double **q_){

    register int j;
    int k;

    for(j=0; j<Nfaces; j++)
        for(k=0; k<Ndof; k++)
            q_[j+Nghost_l][k] = q_[j+Nghost_l][k] + dt_ * resid[j][k];

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
                    + (dt_/6.) * ( resid_temp0[j][k] + 2*resid_temp1[j][k]
                                   + 2*resid_temp2[j][k] + resid[j][k] );

    space_solver->UpdateResid(resid,q_);

    return;
}

void ExplicitTimeSolver::ComputeInitialResid(double **qn_){

    space_solver->UpdateResid(resid,qn_);

    return;
}

void ExplicitTimeSolver::Reset_time_solver(){

    emptyarray(Nfaces,resid);
    emptyarray(Nfaces,q_temp);

    emptyarray(Nfaces,resid_temp0);
    emptyarray(Nfaces,resid_temp1);
    emptyarray(Nfaces,resid_temp2);

    return;
}

