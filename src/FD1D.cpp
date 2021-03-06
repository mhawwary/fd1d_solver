#include "general_tools.h"
#include "global_var.h"
#include"SimData.hpp"
#include"GridData.h"
#include"FDSolver.hpp"
#include"FDSolverAdvec.hpp"
#include"FDSolverAdvecDiffus.hpp"
#include"ExplicitTimeSolver.hpp"
#include"TimeSolver.hpp"


void InitSim(const int& argc, char** argv);
void RunSim();
void PostProcess();
//void Dump_errors_vs_time();
void check_divergence(int& check_div_,int& check_conv_
                      ,const double& threshold, const double& conv_tol
                      , double** Qn_, double& Q_max, double& Q_sum);

void check_convergence(int& check_conv_, double& Q_sum);

void logo();

SimData   simdata;
GridData *meshdata;
FDSolver *fd_solver;
TimeSolver *time_solver;
PadeFilter *filter_solver;

void test_pade_filter();

void test_pade_filter(){

    register int i;
    int nnodes_=25;
    int n_linsys = nnodes_-2;
    double **Q=nullptr;
    double *x=nullptr;

    Q = new double*[nnodes_];
    x = new double[nnodes_];

    double x0=0.0,xf=1.0;
    double dx = (xf-x0)/(nnodes_-1);
    for(i=0; i<nnodes_; i++){
        x[i]   = dx * (i)  + x0 ;  // node 0, element i

        Q[i] = new double[1];
        Q[i][0] = sin(6*PI*x[i]);
    }

    filter_solver = new PadeFilter;
    std::string bound_type = "Periodic";
    filter_solver->setup_filter(simdata.filter_type_,8,nnodes_
                                ,0.40,bound_type);


    printf("\n Unfiltered Q: \n");
    for(i=0; i<nnodes_; i++)
        printf("%1.2f \t %1.5f\n",x[i],Q[i][0]);
    printf("===============================================\n");

    for(i=0; i<25; i++)
        filter_solver->filtered_sol(Q);

    //Q[n_linsys][0] = Q[0][0];

    printf("\n Filtered Q: \n");
    for(i=0; i<nnodes_; i++)
        printf("%1.2f \t %1.10f\n",x[i],Q[i][0]);

    emptyarray(nnodes_,Q);
    emptyarray(x);

    return;
}

int main(int argc, char** argv){

    if (argc < 2) {

        cout << "ERROR: No inputs are specified ... " << endl;

        //test_tridaig(10);
        //test_cyclic_tridiag(10);

        test_pade_filter();

        return 0;
    }

    logo();

    clock_t t_start=clock();

    InitSim(argc, argv);

    RunSim();

    //    if(simdata.wave_form_!=3){
    //        PostProcess();
    //    }else{
    //        printf("\nFinal Iteration number is: %d\n",time_solver->GetIter());
    //        printf("Final time is: %1.5f\n",fd_solver->GetPhyTime());
    //    }

    clock_t t_end=clock();
    cout << "\n\nElapsed Time: " << 1.0*(t_end-t_start)/CLOCKS_PER_SEC
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
        if(simdata.time_solver_type_=="explicit") // Runge-Kutta
            time_solver = new ExplicitTimeSolver;
        time_solver->setupTimeSolver(fd_solver,&simdata);
        fd_solver->dump_timeaccurate_sol();
        fd_solver->dump_timeaccurate_errors();
        simdata.dump_python_inputfile();

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

    /*    int n_iter_print;
//    int local_iter=0;
//    double dt_last_print=0.0;
//    //int check_div_=0,check_conv_=0;
//    //double Q_max=20.0,Q_sum=0.0;

//    double gtime = fd_solver->GetPhyTime();
//    double dt_= fd_solver->GetTimeStep();

//    if(simdata.unsteady_data_print_flag_==0){        // use iter no. to print
//        n_iter_print = simdata.unsteady_data_print_iter_;
//        dt_last_print = dt_;
//        n_iter_print--;
//    }else if(simdata.unsteady_data_print_flag_==1){   // use time point to print
//        n_iter_print= (int) round( simdata.unsteady_data_print_time_/ dt_) ;
//        if((n_iter_print*dt_) > simdata.unsteady_data_print_time_ ){
//            dt_last_print = simdata.unsteady_data_print_time_ - ((n_iter_print-1) * dt_);
//            n_iter_print--;
//        }else if((n_iter_print*dt_) < simdata.unsteady_data_print_time_ ){
//            dt_last_print = simdata.unsteady_data_print_time_ - (n_iter_print*dt_);
//        }
//    }else{
//        FatalError_exit("unsteady data print flag error");
//    }

//    printf("\nnIter to print unsteady data: %d, dt_last: %1.5e"
//           ,n_iter_print, dt_last_print);

//    // First Solve:
//    time_solver->ComputeInitialResid(fd_solver->GetNumSolution());
//    time_solver->SolveOneStep(fd_solver->GetNumSolution());
//    time_solver->space_solver->UpdatePhyTime(dt_);
//    gtime=fd_solver->GetPhyTime();
//    local_iter++;

//    if(time_solver->GetIter()%n_iter_print==0){
//        printf("\nIter No:%d, time: %1.5f\n",time_solver->GetIter(),gtime);
//        fd_solver->dump_timeaccurate_sol();
//    }

    // main solution loop:
//    while ( gtime < (simdata.t_end_-1e-11)  ){

//            time_solver->SolveOneStep(fd_solver->GetNumSolution());
//            time_solver->space_solver->UpdatePhyTime(dt_);
//            gtime=fd_solver->GetPhyTime();
//            local_iter++;

////            check_divergence(check_div_, check_conv_, simdata.div_thresh_, simdata.conv_tol_
////                             ,fd_solver->GetNumSolution(), Q_max, Q_sum);

//            if(local_iter%n_iter_print==0){
//                time_solver->Set_time_step(dt_last_print);
//                time_solver->SolveOneStep(fd_solver->GetNumSolution());
//                time_solver->space_solver->UpdatePhyTime(dt_last_print);
//                gtime=fd_solver->GetPhyTime();
//                printf("\nIter No:%d, time: %1.5f",time_solver->GetIter(),gtime);
//                fd_solver->dump_timeaccurate_sol();
//                time_solver->Set_time_step(dt_);
//                local_iter=0;
//            }

////            if(time_solver->GetIter()%300000==0){
////                cout << "Iter No.:  "<<time_solver->GetIter()
////                     <<"\t Q_max:  "<<Q_max
////                    <<"\t Q_sum:  "<<Q_sum<<endl;
//////                cin.get();
////            }
    }*/

    int n_iter_print;
    int local_iter=0;
    double gtime = fd_solver->GetPhyTime();
    double dt_= fd_solver->GetTimeStep();
    double dt_last_print=0.0;
    double temp_tol=1e-8;

    //======================================================================
    //             Preparing simulation control variables
    //======================================================================
    if(simdata.unsteady_data_print_flag_==0){  // use iter to print
        n_iter_print = simdata.unsteady_data_print_iter_;
        if(n_iter_print<=1)
            FatalError_exit("Warning: iter to print is  <= 1 ");
        dt_last_print = dt_;
        n_iter_print--;

    }else if(simdata.unsteady_data_print_flag_==1){   // use time point to print
        if(simdata.unsteady_data_print_time_ < (dt_ + temp_tol))
            FatalError_exit("Warning:  time to print is less than dt");

        n_iter_print= (int) round( simdata.unsteady_data_print_time_/ dt_) ;
        if((n_iter_print*dt_) > (simdata.unsteady_data_print_time_- temp_tol) ){
            n_iter_print--;
            dt_last_print = simdata.unsteady_data_print_time_ - (n_iter_print * dt_);

        }else if((n_iter_print*dt_) < (simdata.unsteady_data_print_time_+temp_tol) ){
            dt_last_print = simdata.unsteady_data_print_time_ - (n_iter_print*dt_);
        }
        if(n_iter_print<=1)
            FatalError_exit("Warning: iter to print is  <= 1 ");

    // print using the specified iter_print without dt changing except the last one
    }else if(simdata.unsteady_data_print_flag_==2){
        n_iter_print = simdata.unsteady_data_print_iter_;
        dt_last_print=dt_;
    }else{
        FatalError_exit("unsteady data print flag error");
    }

    printf("nIter to print unsteady data: %d, dt_last: %1.5e"
           ,n_iter_print, dt_last_print);
    printf("\n----------------------------------------------------------\n");

    printf("Iter No:%d, time: %f",time_solver->GetIter()
           ,fd_solver->GetPhyTime());

    //===========================
    // Solve First Iteration
    //===========================
    time_solver->ComputeInitialResid(fd_solver->GetNumSolution());
    time_solver->SolveOneStep(fd_solver->GetNumSolution());
    time_solver->space_solver->UpdatePhyTime(dt_);
    gtime=fd_solver->GetPhyTime();
    local_iter++;

    if(n_iter_print==1){
        printf("\nIter No:%d, time: %f",time_solver->GetIter(),gtime);
        fd_solver->dump_timeaccurate_sol();
        fd_solver->dump_timeaccurate_errors();
        //Dump_errors_vs_time(); // for evolution of error vs time
        local_iter=0;
    }

    //======================================================================
    //                        Main Solution Loop
    //======================================================================
    if(simdata.unsteady_data_print_flag_==0
            || simdata.unsteady_data_print_flag_==1){
        while ( fabs(gtime - simdata.t_end_) > (dt_+temp_tol) ){
            time_solver->SolveOneStep(fd_solver->GetNumSolution());
            time_solver->space_solver->UpdatePhyTime(dt_);
            gtime=fd_solver->GetPhyTime();
            local_iter++;

            //check_divergence(check_div_, check_conv_
            // , simdata.div_thresh_, simdata.conv_tol_
            //      ,fd_solver->GetNumSolution(), Q_max, Q_sum

            if(local_iter%n_iter_print==0){
                time_solver->Set_time_step(dt_last_print);
                time_solver->SolveOneStep(fd_solver->GetNumSolution());
                time_solver->space_solver->UpdatePhyTime(dt_last_print);
                gtime=fd_solver->GetPhyTime();
                printf("\nIter No:%d, time: %1.5f",time_solver->GetIter(),gtime);
                fd_solver->dump_timeaccurate_sol();
                fd_solver->dump_timeaccurate_errors();
                //Dump_errors_vs_time(); // for evolution of error vs time
                time_solver->Set_time_step(dt_);
                //time_solver->Reset_iter(time_solver->GetIter()-1);
                local_iter=0;
            }
        }

    }else if(simdata.unsteady_data_print_flag_==2){
        while ( fabs(gtime - simdata.t_end_) > (dt_+temp_tol)){

            time_solver->SolveOneStep(fd_solver->GetNumSolution());
            time_solver->space_solver->UpdatePhyTime(dt_);
            gtime=fd_solver->GetPhyTime();
            local_iter++;

            if(local_iter%n_iter_print==0){
                printf("\nIter No:%d, time: %1.5f",time_solver->GetIter(),gtime);
                fd_solver->dump_timeaccurate_sol();
                fd_solver->dump_timeaccurate_errors();
                //Dump_errors_vs_time(); // for evolution of error vs time
                local_iter=0;
            }
        }

        // Last iteration:
        dt_last_print = fd_solver->GetLastTimeStep();
        time_solver->Set_time_step(dt_last_print);
        time_solver->SolveOneStep(fd_solver->GetNumSolution());
        if(dt_last_print>=temp_tol)
            time_solver->space_solver->UpdatePhyTime(dt_last_print);
        else
            time_solver->space_solver->UpdatePhyTime(dt_);

        gtime=fd_solver->GetPhyTime();
        printf("\nIter No:%d, time: %1.5f",time_solver->GetIter(),gtime);
        fd_solver->dump_timeaccurate_sol();
        fd_solver->dump_timeaccurate_errors();
    }

    return;
}


void PostProcess(){

    double L1_error_=fd_solver->L1_error_nodal_sol();
    double L2_error_=fd_solver->L2_error_nodal_sol();
    double dissip_error = fd_solver->dissipation_error();

    fd_solver->print_cont_vertex_sol();

    fd_solver->dump_errors(L1_error_,L2_error_);

    printf("\nFinal Iteration number is: %d\n",time_solver->GetIter());
    printf("Final time is: %1.2f\n",fd_solver->GetPhyTime());
    printf("nodal_L1_Error: %e , nodal_L2_Error: %e , dissip_Error: %e\n\n"
           ,L1_error_,L2_error_, dissip_error);

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

    cout<<"______________________________________________________________________________________________"<<endl;
    cout<<"                                                                                              "<<endl;
    cout<<"               "<<"  Welcome to the Compact and Finite Difference solver   "<<"               "<<endl;
    cout<<"                   "<<"  for 1D wave and scalar conservation laws  "<<"                       "<<endl;
    cout<<"                                                                                              "<<endl;
    cout<<"         Author:               Mohammad Alhawwary, PhD. Student                               "<<endl;
    cout<<"         Affiliation:   Aerospace Engineering Department, University of Kansas, USA           "<< endl;
    cout<<"__________________________________________03/30/2017__________________________________________"<<endl;

    return;
}


//void Dump_errors_vs_time(){

//    double L1_error_=fd_solver->L1_error_time_nodal_sol();
//    double L2_error_=fd_solver->L2_error_time_nodal_sol();

//    fd_solver->dump_errors_vs_time(L1_error_,L2_error_);

//    printf("\tL1_time_proj: %2.5e\t L2_time_proj: %2.5e",L1_error_,L2_error_);

//    return;
//}






