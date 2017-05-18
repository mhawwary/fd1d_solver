
#include"SimData.hpp"

void SimData::Parse(const std::string &fname){


    GetPot gp_input(fname.c_str());

    Nelem_ = gp_input("Case/num_elements",1);
    x0_ = gp_input("Case/x_begin",0.0);
    xf_ = gp_input("Case/x_end",1.0);
    uniform_ = gp_input("Case/uniform_grid",1);
    refine_level_ = gp_input("Case/refinement_level",0);

    print_freq_=gp_input("Simulation/print_freq",0);
    conv_tol_=gp_input("Simulation/convergence_tolerance",1e-16);
    div_thresh_=gp_input("Simulation/divergence_threshold",20);
    restart_flag = gp_input("Simulation/restart_flag",0);
    restart_iter_ = gp_input("Simulation/restart_iter",0);
    Sim_mode = gp_input("Simulation/mode","normal");

    a_wave_ = gp_input("wave/wave_speed",1);
    wave_form_ = gp_input("wave/wave_form",0);
    Gaussian_exponent_ = gp_input("wave/Gaussian_exponent",-50.0);

    scheme_order_=gp_input("space_solver/order",1);
    upwind_biased_=gp_input("space_solver/upwind_biased",1);

    RK_order_=gp_input("time_solver/explicit/RK_order",0);

    calc_dt_flag = gp_input("time_solver/calculate_dt_flag",1);
    CFL_ = gp_input("time_solver/CFL_no",1e-9);
    dt_ = gp_input("time_solver/dt",1e-5);
    t_init_ = gp_input("time_solver/initial_time",0.0);
    t_end_ = gp_input("time_solver/final_time",1.0);
    Nperiods = gp_input("time_solver/no_of_periods",1.0);
    maxIter_ = gp_input("time_solver/maximum_iteration",1e6);
    end_of_sim_flag_ = gp_input("time_solver/end_of_simulation_flag",1);
}

void SimData::setup_output_directory(){

    /**
    ---------------------------------------------------------
    Creating post processing directory and its subdirectory:
    ---------------------------------------------------------
    */

    struct stat statbuf;

    case_postproc_dir =new char[300];

    char *case_dir=nullptr; case_dir=new char[150];
    char scheme_OA_[20];

    if(scheme_order_==1){
        sprintf(scheme_OA_,"1st");
    }else if(scheme_order_==2){
        sprintf(scheme_OA_,"2nd");
    }else if(scheme_order_==3){
        if(upwind_biased_==1) sprintf(scheme_OA_,"3rd_biased");
        else sprintf(scheme_OA_,"3rd_fupwind");
    }else if(scheme_order_==4){
        sprintf(scheme_OA_,"4th");
    }else{ _notImplemented("Scheme order"); }

    if(Sim_mode=="normal"){
        sprintf(case_dir,"FD%s_RK%d",scheme_OA_,RK_order_);
    }else if(Sim_mode=="test"){
        sprintf(case_dir,"FD%s_RK%d_test",scheme_OA_,RK_order_);
    }else _notImplemented("Simulation mode");

    int test0=0;

    char *current_working_dir=new char[1000];
    getcwd(current_working_dir,1000);

    char results_dir[1500];

    sprintf(results_dir,"%s/Results",current_working_dir);

    if(stat(results_dir, &statbuf) == -1){
          //mkdir("./Results",0777);
        }

    if(!S_ISDIR(statbuf.st_mode)){
        printf("\nResults directory does not exist.....\nCreating Results directory......");
        test0 = mkdir("./Results",0777);
        if(test0==-1) FatalError_exit("Failed to create Results directoy");
    }

    test0 = chdir("./Results");
    if(test0==-1) FatalError_exit("Change directory to ./Results failed");

    test0 = chdir(case_dir);
    if(test0==-1) {
        printf("\nCreating Case_directory.....");
        test0 = mkdir(case_dir,0777);
        if(test0==-1) FatalError_exit("Failed to create case_directory directoy");
        test0 = chdir(case_dir);
        if(test0==-1) FatalError_exit("Change directory to ./Results/case_dir failed");
    }

    sprintf(case_postproc_dir,"./Results/%s/",case_dir);

    mkdir("./nodal",0777);
    mkdir("./errors",0777);
    //mkdir("./BINARY",0777);

    cout<<"\nCurrnet working directory: "<<current_working_dir<<endl;
    cout<<"Post processing directory: "<<case_postproc_dir<<endl;

    test0 = chdir(current_working_dir);
    if(test0==-1) FatalError_exit("Change directory to current_working_dir failed");

    emptyarray(current_working_dir);
    emptyarray(case_dir);

    return;
}

void SimData::dump_python_inputfile(){

    char *fname=nullptr;
    fname = new char[100];

    sprintf(fname,"./input/python_input.in");

    FILE* python_out = fopen(fname,"w");

    fprintf(python_out,"dir:%s\n",case_postproc_dir);
    fprintf(python_out,"exact:%s\n",(char*)"nodal/u_exact");
    fprintf(python_out,"numerical:%s\n",(char*)"nodal/u_num");
    fprintf(python_out,"FDOA:%d\n",scheme_order_);
    fprintf(python_out,"RK:%d\n",RK_order_);
    fprintf(python_out,"Nelem:%d\n",Nelem_);
    fprintf(python_out,"CFL:%1.3f\n",CFL_);
    fprintf(python_out,"T:%1.2f\n",Nperiods);

    fclose(python_out);
    emptyarray(fname);

    return;
}

