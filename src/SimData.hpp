#include"../include/getpot.h"


struct SimData {

    //int restart_flag=0;  //0: start a new simulation; 1: restart a previous simulation
    int scheme_order=0;    // FD Scheme order
    int RK_order=0;        // Runge-Kutta type (0: euler FT, 2: SSPRK22, 3: SSPRK33)
    int print_freq=10;

    double dt = 1e-3;  // dt time step
    double t_init = 0;  // initial time
    double t_end =1e20;  // end time
    double maxIter = 1e10; // maximum number of iterations
    double CFL    = 1.0;   // CFL no.

    int  Nelem = 1;  // no. of elements in the grid

    //std::string case_title;

    //std::string upwind_type_;

    //std::string mesh_fname;   // mesh file name

    //std::string case_postproc_dir;  // postprocessing directory

    void Parse(const std::string &fname);

};
