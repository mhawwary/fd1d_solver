#include "general_tools.h"
#include "global_var.h"
#include"SimData.hpp"
#include"GridData.h"
#include"FDSolver.hpp"
#include"FDSolverAdvec.hpp"
#include"FDSolverAdvecDiffus.hpp"
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

        cout << "ERROR: No inputs are specified ... " << endl;

        //test_tridaig(10);
        //test_cyclic_tridiag(10);

        return 0;
    }

    logo();

    clock_t t_start=clock();

    InitSim(argc, argv);

    RunSim();

    if(simdata.wave_form_!=3)
        PostProcess();

    printf("\nFinal Iteration number is: %d\n",time_solver->GetIter());
    printf("Final time is: %1.5f\n",fd_solver->GetPhyTime());

    clock_t t_end=clock();
    cout << "Elapsed Time: " << 1.0*(t_end-t_start)/CLOCKS_PER_SEC
         << " seconds\n" <<endl;

    emptypointer(meshdata);
    emptypointer(fd_solver);
    emptypointer(time_solver);

    return 0;
}

void InitSim(const int& argc,char** argv){

    if(argc<6){  // Parsing through input file

        simdata.Parse(argv[argc-1]);
        simdata.setup_output_directory();

        if(simdata.wave_form_==3) //Burgers Turbulence
            simdata.prepare_dump_burgers_turb_param();
    }

    meshdata = new GridData;

    if(simdata.restart_flag==0){

        meshdata->set_grid_param(simdata);
        meshdata->generate_grid();

        // Allocating Solvers:
        if(simdata.eqn_set=="Advection")
            fd_solver = new FDSolverAdvec;
        else if(simdata.eqn_set=="Advection_Diffusion")
            fd_solver = new FDSolverAdvecDiffus;
        else
            _notImplemented("Equation set");

        fd_solver->setup_solver(meshdata,simdata);
        fd_solver->InitSol();
        fd_solver->dump_timeaccurate_sol();
        time_solver = new ExplicitTimeSolver;
        time_solver->setupTimeSolver(fd_solver,&simdata);
        simdata.dump_python_inputfile();

        //setup_stencil();
        //init_solution();
        //iter_init = 0;
        //restart_iter_=0;

    } else if(simdata.restart_flag==1){
        printf("\nReading Restart Data\n");
//        BinaryDataReading();
//        generate_grid();
//        setup_stencil();
//        init_sol_fromFile();
//        iter_init = restart_iter_;

    }else {
        FatalError_exit("Wrong restart flag");
    }

    return;
}

void RunSim(){

    int n_iter_print;
    int local_iter=0;
    double dt_last_print=0.0;
    //int check_div_=0,check_conv_=0;
    //double Q_max=20.0,Q_sum=0.0;

    double gtime = fd_solver->GetPhyTime();
    double dt_= fd_solver->GetTimeStep();

    if(simdata.unsteady_data_print_flag_==0){        // use iter no. to print
        n_iter_print = simdata.unsteady_data_print_iter_;
        dt_last_print = dt_;
        n_iter_print--;
    }else if(simdata.unsteady_data_print_flag_==1){   // use time point to print
        n_iter_print= (int) round( simdata.unsteady_data_print_time_/ dt_) ;
        if((n_iter_print*dt_) > simdata.unsteady_data_print_time_ ){
            dt_last_print = simdata.unsteady_data_print_time_ - ((n_iter_print-1) * dt_);
            n_iter_print--;
        }else if((n_iter_print*dt_) < simdata.unsteady_data_print_time_ ){
            dt_last_print = simdata.unsteady_data_print_time_ - (n_iter_print*dt_);
        }
    }else{
        FatalError_exit("unsteady data print flag error");
    }

    printf("\nnIter to print unsteady data: %d, dt_last: %1.5e"
           ,n_iter_print, dt_last_print);

    // First Solve:
    time_solver->ComputeInitialResid(fd_solver->GetNumSolution());
    time_solver->SolveOneStep(fd_solver->GetNumSolution());
    time_solver->space_solver->UpdatePhyTime(dt_);
    gtime=fd_solver->GetPhyTime();
    local_iter++;

    if(time_solver->GetIter()%n_iter_print==0){
        printf("\nIter No:%d, time: %1.5f\n",time_solver->GetIter(),gtime);
        fd_solver->dump_timeaccurate_sol();
    }

    // main solution loop:
    while ( gtime < (simdata.t_end_-1e-11)  ){

            time_solver->SolveOneStep(fd_solver->GetNumSolution());
            time_solver->space_solver->UpdatePhyTime(dt_);
            gtime=fd_solver->GetPhyTime();
            local_iter++;

//            check_divergence(check_div_, check_conv_, simdata.div_thresh_, simdata.conv_tol_
//                             ,fd_solver->GetNumSolution(), Q_max, Q_sum);

            if(local_iter%n_iter_print==0){
                time_solver->Set_time_step(dt_last_print);
                time_solver->SolveOneStep(fd_solver->GetNumSolution());
                time_solver->space_solver->UpdatePhyTime(dt_last_print);
                gtime=fd_solver->GetPhyTime();
                printf("\nIter No:%d, time: %1.5f",time_solver->GetIter(),gtime);
                fd_solver->dump_timeaccurate_sol();
                time_solver->Set_time_step(dt_);
                local_iter=0;
            }

//            if(time_solver->GetIter()%300000==0){
//                cout << "Iter No.:  "<<time_solver->GetIter()
//                     <<"\t Q_max:  "<<Q_max
//                    <<"\t Q_sum:  "<<Q_sum<<endl;
////                cin.get();
//            }
    }

    // Last iteration:

//    time_solver->SolveOneStep(fd_solver->GetNumSolution());
//    time_solver->space_solver->UpdatePhyTime(fd_solver->GetLastTimeStep());
//    gtime=fd_solver->GetPhyTime();

//    printf("Iter No:%d, time: %f\n",time_solver->GetIter(),gtime);
//    fd_solver->dump_timeaccurate_sol();

    return;
}

void PostProcess(){

    double L1_error_=fd_solver->L1_error_nodal_sol();
    double L2_error_=fd_solver->L2_error_nodal_sol();

    fd_solver->print_cont_vertex_sol();

    fd_solver->dump_errors(L1_error_,L2_error_);

    printf("\nFinal Iteration number is: %d\n",time_solver->GetIter());
    printf("Final time is: %1.2f\n",fd_solver->GetPhyTime());
    printf("nodal_L1_Error: %e , nodal_L2_Error: %e\n\n"
           ,L1_error_,L2_error_);

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

    register int i,j;
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









