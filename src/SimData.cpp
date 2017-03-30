#include"SimData.hpp"

void SimData::Parse(const std::string &fname){


    GetPot gp_input(fname.c_str());

    Nelem = gp_input("Case/num_elements",1);

    print_freq=gp_input("Simulation/print_freq",0);

    scheme_order=gp_input("space_solver/order",1);

    RK_order=gp_input("time_solver/explicit/RK_order",0);

    dt = gp_input("time_solver/dt",1e-9);

    t_init = gp_input("time_solver/initial_time",1e-9);

    t_end = gp_input("time_solver/final_time",1e-9);

    maxIter = gp_input("time_solver/maximum_iteration",1e-9);

    CFL = gp_input("time_solver/CFL_no",1e-9);

}
