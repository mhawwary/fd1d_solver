
#include"SimData.hpp"

void SimData::Parse(const std::string &fname){

    GetPot gp_input(fname.c_str());

    // Case parameters:
    //-----------------------------
    case_title = gp_input("Case/title","test");
    Nelem_ = gp_input("Case/num_elements",1);
    x0_ = gp_input("Case/x_begin",0.0);
    xf_ = gp_input("Case/x_end",1.0);
    uniform_ = gp_input("Case/uniform_grid",1);
    refine_level_ = gp_input("Case/refinement_level",0);
    N_exact_ppts = gp_input("Case/N_exact_ploting_pts",100);

    // Simulation parameters:
    //-----------------------------
    unsteady_data_print_flag_=gp_input("Simulation/unsteady_data_print_flag",0);
    unsteady_data_print_iter_=gp_input("Simulation/unsteady_data_print_iter",0);
    unsteady_data_print_time_=gp_input("Simulation/unsteady_data_print_time",0.100);
    conv_tol_=gp_input("Simulation/convergence_tolerance",1e-16);
    div_thresh_=gp_input("Simulation/divergence_threshold",20);
    restart_flag = gp_input("Simulation/restart_flag",0);
    restart_iter_ = gp_input("Simulation/restart_iter",0);
    Sim_mode = gp_input("Simulation/mode","normal");
    case_no_ = gp_input("Simulation/case_no","00");

    // Wave parameters:
    //----------------------------
    a_wave_ = gp_input("wave/wave_speed",1.0);
    wave_form_ = gp_input("wave/wave_form",0);
    // ./trigonometric:
    wave_freq_ = gp_input("wave/trigonometric/wave_frequency",2.0);
    wave_amp_  = gp_input("wave/trigonometric/wave_amplitude",1.0);
    wave_const = gp_input("wave/trigonometric/wave_const",0.0);
    wave_shift = gp_input("wave/trigonometric/wave_freq_shift",0.0);
    // ./Gaussian:
    Gaussian_amp_ = gp_input("wave/Gaussian/Gaussian_amplitude",1.0);
    Gaussian_exponent_ = gp_input("wave/Gaussian/Gaussian_exponent",-50.0);
    // ./Burger_turb:
    if(wave_form_==3){  // Burger's Turbulence
        turb_prob_type_
                = gp_input("wave/Burger_turb/turb_prob_type","Decay_turb_Adams");
        max_wave_no_ = gp_input("wave/Burger_turb/max_wave_no",1024);
        max_energy_wave_no_ = gp_input("wave/Burger_turb/ko",10.0);
        spectrum_restart_flag = gp_input("wave/Burger_turb/spectrum_restart_flag",0);
        velocity_mean_ = gp_input("wave/Burger_turb/velocity_mean",0.0);
    }

    // Space Solver parameters:
    //-----------------------------
    eqn_set = gp_input("space_solver/eqn_set","Advection");
    eqn_type_ = gp_input("space_solver/eqn_type","linear_advec");
    scheme_order_=gp_input("space_solver/order",1);
    // ./advec_eqn:
    upwind_biased_=gp_input("space_solver/advec_eqn/upwind_biased",1);
    // ./heat_eqn:
    thermal_diffus
            = gp_input("space_solver/heat_eqn/thermal_diffusivity",1.0);

    // Time Solver parameters:
    //--------------------------------
    calc_dt_flag = gp_input("time_solver/calculate_dt_flag",1);
    calc_dt_adv_diffus_flag = gp_input("time_solver/calc_dt_adv_diffus_flag",0);
    CFL_ = gp_input("time_solver/CFL_no",1e-9);
    dt_ = gp_input("time_solver/dt",1e-5);
    t_init_ = gp_input("time_solver/initial_time",0.0);
    t_end_ = gp_input("time_solver/final_time",1.0);
    Nperiods = gp_input("time_solver/no_of_periods",1.0);
    maxIter_ = gp_input("time_solver/maximum_iteration",1e6);
    end_of_sim_flag_ = gp_input("time_solver/end_of_simulation_flag",1);
    // ./explicit:
    RK_order_=gp_input("time_solver/explicit/RK_order",0);
}

void SimData::prepare_dump_burgers_turb_param(){

    if(spectrum_restart_flag==0){
        register int i;
        int n_pts_=max_wave_no_;
        double A_=0.;

        k_wave_no_ = new int[n_pts_];
        epsi_phase_ = new double[n_pts_];
        energy_spect_ = new double[n_pts_];

        // Preparing Random number seeds:
        //-----------------------------------
        //srand(time(NULL));
        std::random_device rd;  //Will be used to obtain a seed for the random number engine
        std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
        std::uniform_real_distribution<> dis(0.0
                                             , std::nextafter(0.9999999999999999, DBL_MAX));  // std::nextafter(0.9999999999999999, DBL_MAX)

        // Preparing Dump Spectrum Data:
        //----------------------------------
        char *fname=nullptr;
        fname = new char[400];
        sprintf(fname,"%sspectrum_data_N%d.dat",case_postproc_dir,n_pts_);
        //printf("\nfname_txt_Spect: %s\n",fname);
        FILE* spect_out_ = fopen(fname,"w");

        for(i=0; i<n_pts_; i++){
            epsi_phase_[i] = dis(gen);
            //epsi_phase_[i] = 1.*((double) rand()/(RAND_MAX)) ;
            k_wave_no_[i] = i+1;

            A_ = 2. * pow(max_energy_wave_no_,-5) / (3. * sqrt(PI) ) ;
            energy_spect_[i] = A_ * pow(k_wave_no_[i],4)
                    * exp(-pow((k_wave_no_[i]/max_energy_wave_no_),2)) ;

            fprintf(spect_out_,"%d %2.10e %2.10e\n",k_wave_no_[i]
                    , epsi_phase_[i], energy_spect_[i]);
        }

        fclose(spect_out_);
        emptyarray(fname);

        // Dumping Binary data:
        fname=new char[400];
        sprintf(fname,"%sspectrum_binarydata_N%d",case_postproc_dir,n_pts_);
        //printf("\nfname_bin_Spect: %s\n",fname);
        FILE*  b_spect_out_=fopen(fname,"wb");

        fwrite(&max_wave_no_,sizeof(int),1,b_spect_out_);
        fwrite(k_wave_no_,sizeof(int),n_pts_,b_spect_out_);
        fwrite(epsi_phase_,sizeof(double),n_pts_,b_spect_out_);
        fwrite(energy_spect_,sizeof(double),n_pts_,b_spect_out_);

        fclose(b_spect_out_);
        emptyarray(fname);

    }else if(spectrum_restart_flag==1){    // Reading Binary data
        char *fname=nullptr;
        fname=new char[400];
        sprintf(fname,"%sspectrum_binarydata_N%d",case_postproc_dir,max_wave_no_);
        struct stat statbuf;
        stat(fname, &statbuf);   // check if the binary file exist
        if(!S_ISDIR(statbuf.st_mode)){
            FatalError_exit("Spectrum Binary file does not exist");
        }else{
            FILE*  b_spect_in_=fopen(fname,"rb");
            k_wave_no_ = new int[max_wave_no_];
            epsi_phase_ = new double[max_wave_no_];
            energy_spect_ = new double[max_wave_no_];

            fread(&max_wave_no_,sizeof(int),1,b_spect_in_);
            fread(k_wave_no_,sizeof(int),max_wave_no_,b_spect_in_);
            fread(epsi_phase_,sizeof(double),max_wave_no_,b_spect_in_);
            fread(energy_spect_,sizeof(double),max_wave_no_,b_spect_in_);
            fclose(b_spect_in_);
        }
        emptyarray(fname);
    }else{
        FatalError_exit("spectrum restart flag error");
    }

    return;
}

void SimData::setup_output_directory(){

    /**
    ---------------------------------------------------------
    Creating post processing directory and its subdirectory:
    ---------------------------------------------------------
    */

    struct stat statbuf;

    case_postproc_dir =new char[350];

//    char *case_title=nullptr;
//    case_title=new char[30];
//    if(wave_form_==0) sprintf(case_title,"sine_wave");
//    else if(wave_form_==1) sprintf(case_title,"Gaussian_wave");
//    else if(wave_form_==2) sprintf(case_title,"InViscid_Burgers");
//    else if(wave_form_==3) sprintf(case_title,"Decaying_Burgers_turb");
//    else FatalError_exit("Wrong Wave form and not implemented");

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
    }else if(scheme_order_==6){
        sprintf(scheme_OA_,"6th");
    }else{ _notImplemented("Scheme order"); }

    char *case_dir=nullptr;
    case_dir=new char[70];
    if(Sim_mode=="normal" || Sim_mode=="dt_const" || Sim_mode=="CFL_const"){
        sprintf(case_dir,"FD%s_RK%d",scheme_OA_,RK_order_);
    }else if(Sim_mode=="test"){
        sprintf(case_dir,"FD%s_RK%d_test",scheme_OA_,RK_order_);
    }else _notImplemented("Simulation mode");

    char *main_dir=nullptr;
    main_dir = new char[25];
    if(eqn_set=="Advection"){
        sprintf(main_dir,"Results");
    }else if (eqn_set=="Diffusion"){
        sprintf(main_dir ,"Results_diffus");
    }else if (eqn_set=="Advection_Diffusion"){
        sprintf(main_dir ,"Results_AdvecDiffus");
    }else{
        FatalError_exit("Equation set when specifying output directory");
    }

    char *current_working_dir=new char[1000];
    getcwd(current_working_dir,1000);

    char* results_dir=nullptr;
    results_dir = new char[1050];
    sprintf(results_dir,"%s/%s",current_working_dir,main_dir);

    stat(results_dir, &statbuf);
    int test0=0;
    if(!S_ISDIR(statbuf.st_mode)){
        printf("\nResults directory does not exist.....\nCreating Results directory......");
        //test0 = mkdir("./Results",0777);
        test0 = mkdir(results_dir,0777);
        if(test0==-1) FatalError_exit("Failed to create Results directoy");
    }

    test0 = chdir(results_dir);
    if(test0==-1) FatalError_exit("Change directory to ./Results failed");
    mkdir(case_title.c_str(),0777);
    chdir(case_title.c_str());

    test0 = chdir(case_dir);
    if(test0==-1) {
        printf("\nCreating Case_directory.....");
        test0 = mkdir(case_dir,0777);
        if(test0==-1) FatalError_exit("Failed to create case_directory directoy");
        test0 = chdir(case_dir);
        if(test0==-1) FatalError_exit("Change directory to ./Results/case_dir failed");
    }

    if(wave_form_==3){  // burgers decay turb
        chdir(case_dir);
        char *case_no_t = nullptr;
        case_no_t = new char[100];
        sprintf(case_no_t,"case%s",case_no_.c_str());
        mkdir(case_no_t,0777);
        chdir(case_no_t);
        sprintf(case_no_t,"%s/case%s",case_dir,case_no_.c_str());
        sprintf(case_postproc_dir,"%s/%s/%s/",main_dir,case_title.c_str(),case_no_t);
        emptyarray(case_no_t);
    }else{
        sprintf(case_postproc_dir,"%s/%s/%s/",main_dir,case_title.c_str(),case_dir);
        chdir(case_dir);
    }

    //sprintf(case_postproc_dir,"./%s/%s/%s/",main_dir,case_title,case_dir);

    mkdir("./nodal",0777);
    mkdir("./errors",0777);
    mkdir("./time_data",0777);

    cout<<"\nCurrnet working directory: "<<current_working_dir<<endl;
    cout<<"Post processing directory: "<<case_postproc_dir<<endl;

    test0 = chdir(current_working_dir);
    if(test0==-1) FatalError_exit("Change directory to current_working_dir failed");

    emptyarray(current_working_dir);
    emptyarray(case_dir);
    emptyarray(main_dir);
//    emptyarray(case_title);
    emptyarray(results_dir);

    return;
}

void SimData::dump_python_inputfile(){

    char *fname=nullptr;
    fname = new char[100];

    sprintf(fname,"./input/python_input.in");

    FILE* python_out = fopen(fname,"w");

    fprintf(python_out,"dir:%s\n",case_postproc_dir);
    fprintf(python_out,"errors:%s\n",(char*)"errors/errors");
    fprintf(python_out,"exact:%s\n",(char*)"nodal/u_exact");
    fprintf(python_out,"numerical:%s\n",(char*)"nodal/u_num");
    fprintf(python_out,"unsteady_num:%s\n",(char*)"time_data/u_num");
    fprintf(python_out,"Eqn_set:%s\n",eqn_set.c_str());
    fprintf(python_out,"Eqn_type:%s\n",eqn_type_.c_str());

    if(wave_form_==0) fprintf(python_out,"wave_form:%s\n",(char*)"sine_wave");
    else if(wave_form_==1) fprintf(python_out,"wave_form:%s\n",(char*)"Gaussian_wave");
    else if(wave_form_==2) fprintf(python_out,"wave_form:%s\n",(char*)"InViscid_Burgers");
    else if(wave_form_==3) fprintf(python_out,"wave_form:%s\n",(char*)"Decaying_Burgers_turb");

    fprintf(python_out,"mode:%s\n",Sim_mode.c_str());

    fprintf(python_out,"FDOA:%d\n",scheme_order_);
    fprintf(python_out,"RK:%d\n",RK_order_);
    fprintf(python_out,"Nelem:%d\n",Nelem_);
    fprintf(python_out,"Nexact:%d\n",N_exact_ppts);
    fprintf(python_out,"CFL:%1.3f\n",CFL_);
    fprintf(python_out,"dt:%1.3e\n",dt_);
    fprintf(python_out,"t_end:%1.3e\n",t_end_);
    fprintf(python_out,"T:%1.3f\n",Nperiods);

    fclose(python_out);
    emptyarray(fname);

    return;
}

void SimData::Reset(){

    emptyarray(k_wave_no_);
    emptyarray(epsi_phase_);
    emptyarray(energy_spect_);

    return;
}

