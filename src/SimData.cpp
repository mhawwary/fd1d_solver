#include"SimData.hpp"

void SimData::Parse(const std::string &fname){


    GetPot gp_input(fname.c_str());

    Nelem_ = gp_input("Case/num_elements",1);
    restart_flag = gp_input("Case/restart_flag",0);
    restart_iter_ = gp_input("Case/restart_iter",0);

    print_freq_=gp_input("Simulation/print_freq",0);
    conv_tol_=gp_input("Simulation/convergence_tolerance",1e-16);
    div_thresh_=gp_input("Simulation/divergence_threshold",20);

    scheme_order_=gp_input("space_solver/order",1);
    upwind_biased_=gp_input("space_solver/upwind_biased",1);

    RK_order_=gp_input("time_solver/explicit/RK_order",0);

    dt_ = gp_input("time_solver/dt",1e-5);

    t_init_ = gp_input("time_solver/initial_time",0.0);

    maxtime_ = gp_input("time_solver/final_time",1.0);

    maxIter_ = gp_input("time_solver/maximum_iteration",1e6);

    max_iter_flag_ = gp_input("time_solver/max_iter_flag",1);

    CFL_ = gp_input("time_solver/CFL_no",1e-9);

    a_wave_ = gp_input("wave/wave_speed",1);

}

