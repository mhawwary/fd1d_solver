#include "general_tools.h"
#include "global_var.h"
//#include "space_solve.h"
//#include "time_solve.h"
#include "post_process.h"
#include "solver_func.h"
#include "grid.h"


void read_input(char** argv_);
void setsizes();
void InitSim(char** argv);
void RunSim();
void PostProcess();


int main(int argc, char** argv){

    if (argc < 2) {

        cout << "ERROR: No inputs are specified ... " << endl; return(0);
    }

    InitSim(argv);

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

    Qn = new double[ Nfaces];

    Qtemp = new double [ Nghost_l + Nfaces + Nghost_r ];

    return;

}

void InitSim(char** argv){

    read_input(argv);

    setsizes();

    generate_grid();

    setup_stencil();

    register int i;

    for( i=0; i<Nfaces; i++ ){

        Qinit[i] = sin(2*PI*X[i]);

        Qn[i] = Qinit[i];

        //_print(i,Qn[i]);

        //update_ghost_sol();
    }

    //initial_ghost_sol_dump();
    initial_solution_dumping();

    return;
}

void RunSim(){

    int n=0;

    gtime=0;

    while ( gtime <=fabs( max_time - pow(10,-10)) ){

        gtime += dt;  n++;

        ComputeOneStep();

        //initial_ghost_sol_dump();

        //intermediate_solution_dump(n, gtime);
    }

    cout << "\nNo. of Iterations:  "<<n<<endl;
    cout << "\nActual time:  "<<max_time<<endl;

    return;
}


void PostProcess(){

    double L1_norm=0.0;
    double L2_norm=0.0;

    final_solution_dump();

    //ComputeError(L1_norm, L2_norm);

    //serror_dumping(L1_norm, L2_norm);

    return;
}














