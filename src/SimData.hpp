#ifndef SIMDATA_H
#define SIMDATA_H

#include"../include/getpot.h"
#include"general_tools.h"
#include"global_var.h"


struct SimData {

    std::string eqn_set;

    int  Nelem_ = 1;  // no. of elements in the grid
    double x0_=0.0;
    double xf_=1.0;
    int uniform_=1;  // flag for uniform/nonunifrom grids, 1: uniform
    int refine_level_=0; // 0: no refinement
    int N_exact_ppts = 100;

    int print_freq_=10;
    double conv_tol_=1e-15;
    double div_thresh_ = 20;
    int restart_flag =0;  //0: start a new simulation; 1: restart a previous simulation
    unsigned int restart_iter_=100;
    std::string Sim_mode;
    std::string case_no_;  // case no for burgers decay turb

    double a_wave_=2;
    int wave_form_ = 0;  // 0: sine wave, 1: Gaussian wave
    double wave_freq_=2.0;  // wave_freq_ x PI
    double Gaussian_exponent_ = -40; // u(x) = exp(-38.6 *x^2)
    double thermal_diffus=1.0;

    int scheme_order_=0;    // FD Scheme order
    int RK_order_=0;        // Runge-Kutta type (0: euler FT, 2: SSPRK22, 3: SSPRK33)
    int upwind_biased_=1;

    int calc_dt_flag=1; // 1: specify CFL and calc dt, 0: specify dt and calc CFL
    double CFL_    = 1.0;   // CFL no.
    double dt_ = 1e-3;  // dt time step
    double t_init_ = 0;  // initial time
    double t_end_ =1e20;  // end time
    unsigned int maxIter_ = 1e6; // maximum number of iterations
    double Nperiods = 1.0; // no. of periods for simulation
    int end_of_sim_flag_=0;  // 1: use max_iteration as a stopping criteria if not converged or diverged

    // Burger's Tubulence Parameters:
    std::string turb_prob_type_;
    int max_wave_no_ = 1024;
    double max_energy_wave_no_ = 10.0;
    int* k_wave_no_ =nullptr;
    double* epsi_phase_=nullptr;
    double* energy_spect_=nullptr;
    int spectrum_restart_flag=0;
    double data_print_time_=0.01;

    char* case_postproc_dir=nullptr;  // postprocessing directory

    void Parse(const std::string &fname);

    void setup_output_directory();

    void dump_python_inputfile();

    void prepare_dump_burgers_turb_param();

    void Reset();
};

#endif
