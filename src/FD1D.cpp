#include "general_tools.h"
#include "global_var.h"
#include"SimData.hpp"
#include"GridData.h"
#include"FDSolver.hpp"
#include"ExplicitTimeSolver.hpp"

void InitSim(const int& argc, char** argv);
void RunSim();
void PostProcess();
void check_divergence(int& check_div_,int& check_conv_
                      ,const double& threshold, const double& conv_tol
                      , double** Qn_, double& Q_max, double& Q_sum);

void check_convergence(int& check_conv_, double& Q_sum);

void logo();

SimData   simdata;
GridData *meshdata;
FDSolver *fd_solver;
ExplicitTimeSolver *time_solver;


int main(int argc, char** argv){

    if (argc < 2) {

        cout << "ERROR: No inputs are specified ... " << endl; return(0);
    }

    logo();

    InitSim(argc, argv);

    RunSim();

    PostProcess();

    emptypointer(meshdata);
    emptypointer(fd_solver);
    emptypointer(time_solver);

    return 0;
}

void InitSim(const int& argc,char** argv){

    if(argc<6){  // Parsing through input file

        simdata.Parse(argv[argc-1]);
        simdata.setup_output_directory();
        simdata.dump_python_inputfile();
    }

    meshdata = new GridData;

    if(simdata.restart_flag==0){

        meshdata->set_grid_param(simdata);
        meshdata->generate_grid();

        fd_solver= new FDSolver;

        fd_solver->setup_solver(meshdata,simdata);

        fd_solver->InitSol();

        time_solver = new ExplicitTimeSolver;

        //setup_stencil();

        //init_solution();

        //iter_init = 0;

        //restart_iter_=0;

    } else if(simdata.restart_flag==1){

//        BinaryDataReading();

//        generate_grid();

//        setup_stencil();

//        init_sol_fromFile();

//        iter_init = restart_iter_;

    }else {
        FatalError("Wrong restart flag");
    }

    return;
}

void RunSim(){

//    int check_div_=0,check_conv_=0;
//    double Q_max=20.0,Q_sum=0.0;

    time_solver->setupTimeSolver(fd_solver,&simdata);

    double gtime = fd_solver->GetPhyTime();

    double dt_= fd_solver->GetTimeStep();

    time_solver->ComputeInitialResid(fd_solver->GetNumSolution());

    time_solver->SolveOneStep(fd_solver->GetNumSolution());

    time_solver->space_solver->UpdatePhyTime(dt_);

    gtime=fd_solver->GetPhyTime();

    while ( gtime < (simdata.t_end_-0.5*dt_) ){

            time_solver->SolveOneStep(fd_solver->GetNumSolution());

            time_solver->space_solver->UpdatePhyTime(dt_);

            gtime=fd_solver->GetPhyTime();

//            check_divergence(check_div_, check_conv_, simdata.div_thresh_, simdata.conv_tol_
//                             ,fd_solver->GetNumSolution(), Q_max, Q_sum);

//            if(time_solver->GetIter()%300000==0){
//                cout << "Iter No.:  "<<time_solver->GetIter()
//                     <<"\t Q_max:  "<<Q_max
//                    <<"\t Q_sum:  "<<Q_sum<<endl;
//            }
    }

    return;
}

void PostProcess(){

    double L2_norm_error=fd_solver->ComputeSolNodalError();

    fd_solver->print_cont_vertex_sol();

    fd_solver->dump_errors(L2_norm_error);

    printf("\nFinal Iteration number is: %d\n",time_solver->GetIter());
    printf("Final time is: %1.2f\n",fd_solver->GetPhyTime());
    printf("nodal_L2_Error: %e\n\n",L2_norm_error);

    return;
}


void check_divergence(int& check_div_, int& check_conv_,const double& threshold
                      , const double& conv_thresh,double** Qn_
                      , double& Q_max, double& Q_sum){

    register int i;

    int Nghost_l = fd_solver->GetNghost_l();
    int Nfaces = fd_solver->GetNfaces();

    Q_max = fabs(Qn_[Nghost_l][0]);
    Q_sum = fabs(Qn_[Nghost_l][0]);

    for(i=Nghost_l+1; i<Nfaces; i++){

        if(fabs(Qn_[i][0]) > Q_max) Q_max = fabs(Qn_[i][0]);

        Q_sum+= fabs(Qn_[i][0]);
    }

    if(Q_max >= threshold) check_div_=1;
    if(Q_max <= conv_thresh) check_conv_=1;

    return;
}

void check_convergence(int& check_conv_, double& Q_sum){

//    register int i,j;

//    Q_sum = 0.0;

//    for(i=Nghost_l; i<Nfaces; i++){

//        Q_sum+= fabs(Qn[i]);
//    }

//    if(Q_sum <= 1e-10) check_conv_=1;


    return;
}


void logo(){

    cout<<"_________________________________________________________________________________________"<<endl;
    cout<<"                                                                                         "<<endl;
    cout<<"                "<<"  Welcome to the Finite Difference solver   "<<"                  "<<endl;
    cout<<"                "<<"  for 1D wave and scalar conservation laws  "<<"                  "<<endl;
    cout<<"                                                                                         "<<endl;
    cout<<"         Author:               Mohammad Alhawwary, PhD. Student                          "<<endl;
    cout<<"    Affiliation:   Aerospace Engineering Department, University of Kansas, USA           "<< endl;
    cout<<"_______________________________________03/30/2017________________________________________"<<endl;

    return;
}









