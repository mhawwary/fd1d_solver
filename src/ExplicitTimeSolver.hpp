#ifndef EXPLICITTIMESOLVER_H
#define EXPLICITTIMESOLVER_H

#include"TimeSolver.hpp"
#include"FDSolver.hpp"
#include"SimData.hpp"

class ExplicitTimeSolver:public TimeSolver{

public:
    ExplicitTimeSolver(void);
    virtual ~ExplicitTimeSolver(void);
    virtual void setupTimeSolver(FDSolver* fd_solver_, SimData* simdata_);
    virtual void SolveOneStep(double **qn_);
    virtual void ComputeInitialResid(double **qn_);
    virtual void Reset_time_solver();

protected:
    void FwdEuler(double **q_);
    void SSPRK22(double **q_);
    void SSPRK33(double **q_);
    void classicRK4(double **q_);

    void leapfrog(double **q_); // standard leapfrog scheme
    void upwind_leapfrog(double **q_); // upwind leapfrog scheme, Isreles86 , Thomas&Roe93, Kim&Roe1999

    void CopyOldSol(double **q_t_, double **qn_);
    void CopyOldestSol(double **q_t_, double **qn_);
    void CopyOldSolwithGhost(double **q_t_, double **qn_);
    void CopyOldResid(double **resid_t_, double **old_resid_);

protected:
    double **q_temp=nullptr;
    double **qnm1_temp=nullptr;
    double **resid_temp0=nullptr;
    double **resid_temp1=nullptr;
    double **resid_temp2=nullptr;

};

#endif
