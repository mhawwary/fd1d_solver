#include "general_tools.h"
//#include "SimData.hpp"
//#include "GridData.hpp"
#include "global_var.h"
#include"space_solve.h"
#include"time_solve.h"
//#include"invflux.h"
#include"post_process.h"


void read_input(char** argv_);
void setsizes();
void InitSim();
void RunSim();
void PostProcess();


int main(int argc, char** argv){

    if (argc < 2) {

        cout << "ERROR: No inputs are specified ... " << endl; return(0);
    }

    read_input(argv);

    generate_grid();

    InitSim();

    RunSim();


//    std::string input_fname;     // input file name

//    input_fname = argv[argc-1];  // input file name with directory

//    SimData  simdata_;

//    GridData meshdata_;

//    simdata_.Parse(input_fname);  // Setup and parse input parameters

//    InitSim(simdata_, meshdata_);           // Preprocessing steps

//    RunSim(simdata_, meshdata_);            // Main Solution Loop

//    PostProcess(simdata_, meshdata_);       // Dumping Simulation Post Processing data


    return 0;
}


void read_input(char** argv){

    Nelem = atof(argv[1]);

    // Sample prompt input:  Nelem Max_time cflno./dt AA
    Nelem = atof(argv[1]);      // Number of elements
    max_time = atof(argv[2]);   // Max time for end of simulation
    CFL = atof(argv[3]);      // CFL number
    //dt = atof(argv[3]);       // time step
    a_wave    = atof(argv[4]);      // Wave speed

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

    return;
}

void setsizes(){

    register int i;

    X     = new double[Nfaces];

    h_j = new double[Nelem];
    local_cfl = new double[Nelem];

    Qinit  = new double[Nfaces];
    Qn  = new double[Nfaces];
    Qex = new double[Nfaces];

    Resid = new double [Nfaces];

    stencil_index = new int[scheme_order+1];
    FD_coeff = new int[scheme_order+1];

    return;

}

void InitSim(){

    register int i,j;

    for( i=0; i<Nfaces; i++ )
        Qinit[i] = sin(PI*X[i]);

    //initial_solution_dumping();

    return;
}

void RunSim(){

    setup_stencil();

    int n=0;

    gtime=0;

    while ( gtime <=fabs( max_time - pow(10,-10)) ){

        gtime += dt;  n++;

        ComputeOneStep();
    }


    return;
}


