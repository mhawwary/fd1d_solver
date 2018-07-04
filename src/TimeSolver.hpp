#ifndef TIMESOLVER_H
#define TIMESOLVER_H

#include"FDSolver.hpp"
#include"SimData.hpp"

class TimeSolver{

public:
    TimeSolver(void){}
    virtual ~TimeSolver(void){}
    virtual void setupTimeSolver(FDSolver* fd_solver_, SimData* simdata_)=0;
    virtual void SolveOneStep(double **qn_)=0;
    virtual void ComputeInitialResid(double **qn_)=0;

    int GetIter(){
        return IterNo;
    }

    void Set_time_step(const double& dt_set){
        dt_ = dt_set;
        return;
    }

    void Reset_iter(const double& iter_reset){
        IterNo = iter_reset;
        return;
    }

protected:
    void UpdateIter(){

        IterNo++;

        return;
    }

    virtual void Reset_time_solver()=0;

public:
    FDSolver *space_solver=nullptr;
    SimData  *simdata=nullptr;

protected:
    std::string time_scheme_type;
    double **resid=nullptr;

    int n_time_levels=1;  // number of time levels
    int Nfaces=1;

    int Nghost_l=1;

    int Ndof=1;

    int IterNo=0;

    double dt_=0.0;

};

#endif
