#include "general_tools.h"
#include "global_var.h"
#include"SimData.hpp"
#include"GridData.h"
#include"FDSolver.hpp"
#include"ExplicitTimeSolver.hpp"

//void read_input(char** argv_);
//void setsizes();
void InitSim(const int& argc, char** argv);
void RunSim();
void PostProcess();
void check_divergence(int& check_div_,int& check_conv_
                      ,const double& threshold, const double& conv_tol
                      , double& Q_max, double& Q_sum);

void check_convergence(int& check_conv_, double& Q_sum);

void logo();

GridData *meshdata;
SimData   simdata;
FDSolver *fd_solver;
ExplicitTimeSolver *time_solver;


int main(int argc, char** argv){

    if (argc < 2) {

        cout << "ERROR: No inputs are specified ... " << endl; return(0);
    }

    logo();

    InitSim(argc, argv);

    RunSim();

    //PostProcess();

    emptypointer(meshdata);
    emptypointer(fd_solver);
    emptypointer(time_solver);

    return 0;
}

//void read_input(char** argv){

//    // Sample prompt input:  Nelem Max_time cflno./dt AA
//    Nelem = atof(argv[1]);      // Number of elements
//    max_time = atof(argv[2]);   // Max time for end of simulation
//    CFL = atof(argv[3]);      // CFL number
//    //dt = atof(argv[3]);       // time step
//    a_wave    = atof(argv[4]);      // Wave speed
//    scheme_order = atof(argv[5]); // scheme order
//    RK_order = atof(argv[6]); // RK order

//    Nfaces = Nelem + 1;
//    Ndof = Nelem;
//    dx=(xf-x0)/Nelem;
//    dt=dx*CFL/a_wave;
//    //sigma = AA*dt/dx;

//    // Screen Output of input and simulation parameters:
//    cout << "CFL no.:  "<<CFL<<"\tWave Speed:  "<<a_wave<<endl;
//    cout << "dt:  "<<dt<<"\t"<< "dx:  "<<dx<<endl;
//    cout << "required no. of time steps: "<<max_time/dt<<endl;
//    cout << "Number of Elements:  "<<Nelem<<endl;
//    cout << "Scheme Order: "<<scheme_order<<endl;
//    cout << "RK_order:  "<< RK_order << endl;

//    return;
//}

//void setsizes(){

//    X     = new double[Nfaces];

//    h_j = new double[Nelem];
//    local_cfl = new double[Nelem];

//    Qinit  = new double[Nfaces];
//    Qex = new double[Nfaces];

//    Resid = new double [Nfaces];

//    stencil_index = new int[scheme_order+1];
//    FD_coeff = new double [scheme_order+1];

//    if(scheme_order==1){
//        Nghost_l = 1;
//        Nghost_r=0;

//    }else if(scheme_order==2){
//        Nghost_l=1;
//        Nghost_r=1;

//    }else if(scheme_order==3){
//        Nghost_l=2;
//        Nghost_r=1;

//    }else if(scheme_order==4){
//        Nghost_l=2;
//        Nghost_r=2;
//    }

//    Netot = Nghost_l + Nfaces + Nghost_r ;
//    Qn = new double[ Netot];

//    return;

//}

void InitSim(const int& argc,char** argv){

    if(argc<6){  // Parsing through input file

        simdata.Parse(argv[argc-1]);
        simdata.setup_output_directory();
    }

    meshdata = new GridData;

    if(simdata.restart_flag==0){

        meshdata->set_grid_param(simdata);
        meshdata->generate_grid();

        fd_solver= new FDSolver;

        fd_solver->setup_solver(meshdata,simdata);

        fd_solver->InitSol();

        fd_solver->print_cont_vertex_sol();
        double L2_error =fd_solver->ComputeSolNodalError();
        fd_solver->dump_errors(L2_error);

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

    time_solver->setupTimeSolver(fd_solver,&simdata);

    double gtime = fd_solver->GetPhyTime();

    double dt_= fd_solver->GetTimeStep();

    //printf("Address of Qn before sending to timesolver: %d\n",fd_solver->GetNumSolution());

    time_solver->ComputeInitialResid(fd_solver->GetNumSolution());

    time_solver->SolveOneStep(fd_solver->GetNumSolution());

    time_solver->space_solver->UpdatePhyTime(dt_);

    gtime=fd_solver->GetPhyTime();

    while ( gtime < (simdata.t_end_-0.5*dt_) ){

            time_solver->SolveOneStep(fd_solver->GetNumSolution());

            time_solver->space_solver->UpdatePhyTime(dt_);

            gtime=fd_solver->GetPhyTime();
    }

//    unsigned int n=iter_init;

//    int check_div_=0,check_conv_=0;

//    gtime=t_init;

//    double Q_max=20.0,Q_sum=0.0;

//    if(max_iter_flag==0){

//        while ( gtime <=fabs( max_time - pow(10,-10)) ){

//            gtime += dt;  n++;

//            ComputeOneStep();

//            check_divergence(check_div_, check_conv_, div_thresh_, conv_tol_, Q_max, Q_sum);

//            //check_convergence(check_conv_, Q_sum);

//            if((n%1000000)==0) {

//                cout << "Iter No.:  "<<n
//                     <<"\t Q_max:  "<<Q_max
//                    <<"\t Q_sum:  "<<Q_sum<<endl;

//                BinaryDataWriting(n);
//            }

//            if(check_div_ ==1) { cout <<"++++++ Diverged +++++++\tQ_max:  "<<Q_max<<"\n\n";  break; }
//            if(check_conv_==1) { cout <<"++++++ Converged +++++++\tQ_max:  "<<Q_max<<"\n\n"; break; }

//        }

//    }else {

//        while ( n < ( max_iter_ + restart_iter_)  ){

//            gtime += dt;  n++;

//            ComputeOneStep();

//            check_divergence(check_div_, check_conv_, div_thresh_, conv_tol_, Q_max, Q_sum);

//            //check_convergence(check_conv_, Q_sum);

//            if((n%1000000)==0) {

//                cout << "Iter No.:  "<<n
//                     <<"\t Q_max:  "<<Q_max
//                    <<"\t Q_sum:  "<<Q_sum<<endl;

//                BinaryDataWriting(n);
//            }

//            if(check_div_ ==1) { cout <<"++++++ Diverged +++++++\tQ_max:  "<<Q_max<<"\n\n";  break; }
//            if(check_conv_==1) { cout <<"++++++ Converged +++++++\tQ_max:  "<<Q_max<<"\n\n"; break; }

//        }

//    }



//    cout << "No. of Iterations:  "<<n<<endl;
//    cout << "Actual time:  "<<max_time<<endl;
//    cout <<"\n===============================================\n";

    return;
}

void PostProcess(){

//    double L1_norm=0.0;
//    double L2_norm=0.0;

//    Compute_exact_shifted_sol();

//    final_solution_dump();

//    ComputeError(L1_norm, L2_norm);

//    //error_dumping(L1_norm, L2_norm);

    return;
}


void check_divergence(int& check_div_, int& check_conv_,const double& threshold, const double& conv_thresh, double& Q_max, double& Q_sum){

//    register int i;

//    Q_max = fabs(Qn[Nghost_l]);
//    Q_sum = fabs(Qn[Nghost_l]);

//    for(i=Nghost_l+1; i<Nfaces; i++){

//        if(fabs(Qn[i]) > Q_max) Q_max = fabs(Qn[i]);

//        Q_sum+= fabs(Qn[i]);
//    }

//    if(Q_max >= threshold) check_div_=1;
//    if(Q_max <= conv_thresh) check_conv_=1;

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









