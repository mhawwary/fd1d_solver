#include "general_tools.h"
#include "global_var.h"
#include "post_process.h"
#include "solver_func.h"
#include "grid.h"
#include"SimData.hpp"

void update_globalVar(SimData& simdata_);
void read_input(char** argv_);
void setsizes();
void InitSim(const int& argc, char** argv);
void RunSim();
void PostProcess();
void check_divergence(int& check_div_,int& check_conv_
                      ,const double& threshold, const double& conv_tol
                      , double& Q_max, double& Q_sum);

void check_convergence(int& check_conv_, double& Q_sum);

void setup_output_directory();

int main(int argc, char** argv){

    if (argc < 2) {

        cout << "ERROR: No inputs are specified ... " << endl; return(0);
    }

    setup_output_directory();

    InitSim(argc, argv);

    RunSim();

    PostProcess();

    return 0;
}

void read_input(char** argv){

    // Sample prompt input:  Nelem Max_time cflno./dt AA
    Nelem = atof(argv[1]);      // Number of elements
    max_time = atof(argv[2]);   // Max time for end of simulation
    CFL = atof(argv[3]);      // CFL number
    //dt = atof(argv[3]);       // time step
    a_wave    = atof(argv[4]);      // Wave speed
    scheme_order = atof(argv[5]); // scheme order
    RK_order = atof(argv[6]); // RK order

    Nfaces = Nelem + 1;
    Ndof = Nelem;
    dx=(xf-x0)/Nelem;
    dt=dx*CFL/a_wave;
    //sigma = AA*dt/dx;

    // Screen Output of input and simulation parameters:
    cout << "CFL no.:  "<<CFL<<"\tWave Speed:  "<<a_wave<<endl;
    cout << "dt:  "<<dt<<"\t"<< "dx:  "<<dx<<endl;
    cout << "required no. of time steps: "<<max_time/dt<<endl;
    cout << "Number of Elements:  "<<Nelem<<endl;
    cout << "Scheme Order: "<<scheme_order<<endl;
    cout << "RK_order:  "<< RK_order << endl;

    return;
}

void setsizes(){

    X     = new double[Nfaces];

    h_j = new double[Nelem];
    local_cfl = new double[Nelem];

    Qinit  = new double[Nfaces];
    Qex = new double[Nfaces];

    Resid = new double [Nfaces];

    stencil_index = new int[scheme_order+1];
    FD_coeff = new double [scheme_order+1];

    if(scheme_order==1){
        Nghost_l = 1;
        Nghost_r=0;

    }else if(scheme_order==2){
        Nghost_l=1;
        Nghost_r=1;

    }else if(scheme_order==3){
        Nghost_l=2;
        Nghost_r=1;

    }else if(scheme_order==4){
        Nghost_l=2;
        Nghost_r=2;
    }

    Netot = Nghost_l + Nfaces + Nghost_r ;
    Qn = new double[ Netot];

    return;

}

void InitSim(const int& argc,char** argv){

    SimData simdata_;

    if(argc<6){  // Parsing through input file

        simdata_.Parse(argv[argc-1]);
        update_globalVar(simdata_);

    }else {
        read_input(argv);
    }

    setsizes();

    if(simdata_.restart_flag==0){

        generate_grid();

        setup_stencil();

        init_solution();

        iter_init = 0;

        restart_iter_=0;

    } else if(simdata_.restart_flag==1){

        BinaryDataReading();

        generate_grid();

        setup_stencil();

        init_sol_fromFile();

        iter_init = restart_iter_;

    }else {
        FatalError("Wrong restart flag");
    }

    return;
}

void update_globalVar(SimData& simdata_){

    Nelem = simdata_.Nelem_;

    CFL= simdata_.CFL_;

    dt = simdata_.dt_;

    max_time = simdata_.maxtime_;

    RK_order = simdata_.RK_order_;

    scheme_order = simdata_.scheme_order_;

    upwind_biased = simdata_.upwind_biased_;

    a_wave = simdata_.a_wave_;

    Nfaces = Nelem + 1;
    Ndof = Nelem;
    dx=(xf-x0)/Nelem;
    dt=dx*CFL/a_wave;

    conv_tol_ = simdata_.conv_tol_;
    div_thresh_ = simdata_.div_thresh_;
    restart_iter_ = simdata_.restart_iter_;
    max_iter_ = simdata_.maxIter_;
    max_iter_flag = simdata_.max_iter_flag_;


    // Screen Output of input and simulation parameters:
    cout <<"\n===============================================\n";
    cout << "CFL no.:  "<<CFL<<"\tWave Speed:  "<<a_wave<<endl;
    cout << "dt:  "<<dt<<"\t"<< "dx:  "<<dx<<endl;
    cout << "required no. of time steps: "<<max_time/dt<<endl;
    cout << "Number of Elements:  "<<Nelem<<endl;
    cout << "Scheme Order: "<<scheme_order<<endl;
    cout << "RK_order:  "<< RK_order << endl <<"\n";

    return;
}

void RunSim(){

    unsigned int n=iter_init;

    int check_div_=0,check_conv_=0;

    gtime=t_init;

    double Q_max=20.0,Q_sum=0.0;

    if(max_iter_flag==0){

        while ( gtime <=fabs( max_time - pow(10,-10)) ){

            gtime += dt;  n++;

            ComputeOneStep();

            check_divergence(check_div_, check_conv_, div_thresh_, conv_tol_, Q_max, Q_sum);

            //check_convergence(check_conv_, Q_sum);

            if((n%1000000)==0) {

                cout << "Iter No.:  "<<n
                     <<"\t Q_max:  "<<Q_max
                    <<"\t Q_sum:  "<<Q_sum<<endl;

                BinaryDataWriting(n);
            }

            if(check_div_ ==1) { cout <<"++++++ Diverged +++++++\tQ_max:  "<<Q_max<<"\n\n";  break; }
            if(check_conv_==1) { cout <<"++++++ Converged +++++++\tQ_max:  "<<Q_max<<"\n\n"; break; }

        }

    }else {

        while ( n < ( max_iter_ + restart_iter_)  ){

            gtime += dt;  n++;

            ComputeOneStep();

            check_divergence(check_div_, check_conv_, div_thresh_, conv_tol_, Q_max, Q_sum);

            //check_convergence(check_conv_, Q_sum);

            if((n%1000000)==0) {

                cout << "Iter No.:  "<<n
                     <<"\t Q_max:  "<<Q_max
                    <<"\t Q_sum:  "<<Q_sum<<endl;

                BinaryDataWriting(n);
            }

            if(check_div_ ==1) { cout <<"++++++ Diverged +++++++\tQ_max:  "<<Q_max<<"\n\n";  break; }
            if(check_conv_==1) { cout <<"++++++ Converged +++++++\tQ_max:  "<<Q_max<<"\n\n"; break; }

        }

    }



    cout << "No. of Iterations:  "<<n<<endl;
    cout << "Actual time:  "<<max_time<<endl;
    cout <<"\n===============================================\n";

    return;
}

void PostProcess(){

    double L1_norm=0.0;
    double L2_norm=0.0;

    Compute_exact_shifted_sol();

    final_solution_dump();

    ComputeError(L1_norm, L2_norm);

    _print(L1_norm,L2_norm);
    cout <<"                         \n";

    //error_dumping(L1_norm, L2_norm);

    return;
}


void check_divergence(int& check_div_, int& check_conv_,const double& threshold, const double& conv_thresh, double& Q_max, double& Q_sum){

    register int i;

    Q_max = fabs(Qn[Nghost_l]);
    Q_sum = fabs(Qn[Nghost_l]);

    for(i=Nghost_l+1; i<Nfaces; i++){

        if(fabs(Qn[i]) > Q_max) Q_max = fabs(Qn[i]);

        Q_sum+= fabs(Qn[i]);
    }

    if(Q_max >= threshold) check_div_=1;
    if(Q_max <= conv_thresh) check_conv_=1;

    return;
}


void check_convergence(int& check_conv_, double& Q_sum){

    register int i,j;

    Q_sum = 0.0;

    for(i=Nghost_l; i<Nfaces; i++){

        Q_sum+= fabs(Qn[i]);
    }

    if(Q_sum <= 1e-10) check_conv_=1;


    return;
}

void setup_output_directory(){

    /**
    ---------------------------------------------------------
    Creating post processing directory and its subdirectory:
    ---------------------------------------------------------
    */
    allocator<char> allchar; // default allocator for char

    char *current_working_dir=allchar.allocate(1000);

    getcwd(current_working_dir,1000);

    chdir("./output");

    mkdir("./BINARY",0777);

    chdir(current_working_dir);

    cout<<"--> Currnet working directory: "<<current_working_dir<<endl;

    return;
}











