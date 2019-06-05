///      main.cpp -- Aggregation Orientation model
///      Date: 2012
///      Mikael Mohtaschemi

#include <iostream>
#include <cstdlib>
#include "Combo_class.hpp"

#include "cvode/cvode.h"
#include "cvode/cvode_dense.h"
#include "cvode/cvode_lapack.h"
#include "nvector/nvector_serial.h"
#include "sundials/sundials_dense.h"
#include "sundials/sundials_types.h"
#include "sundials/sundials_math.h"
#include "sundials/sundials_spgmr.h"
#include "cvodes/cvodes_spgmr.h"

using namespace std;
int wfunc (realtype t, N_Vector y, N_Vector ydot, void* userdata);
///THE MAIN FUNCTION
int main(int argc, char **argv)
{
    ///initialize stuff
    double t0 = 0, dt = 1e-6, time = 0, max_time = 100;
    int s1, s2;
    ///diffusion coefficient, aspect ratio, shear rate, beta=v_d/v_c
    double diffusion_coeff, aspect_ratio, shear_rate, beta=1.0;
    int theta_grid, phi_grid, classes=1;
    int what_todo;
    cout << "###########################################################\n";
    cout << "############ AGGREGATION ORIENTATION MODEL ################\n";
    cout << "###########################################################\n";
    cout << "#INPUT: executable phi_grid theta_grid diffusion_coeff ";
    cout << "shear_rate  aspect_ratio dt max_time s1 s2 classes what_todo visco_0 Np_max beta shear_rate_max" << endl;
    cout << "what_todo::\n-1 = selfcheck\n 0=new run\n 1=continue from last state" << endl; 
    double visco0 = 1.0; //Some default values for this
    double Np_max = 1.0; //and this
    double shear_rate_max = 100.0;
    if(argc >= 16){
		shear_rate_max = stof(argv[15]);
	}    
    if(argc >=15){
        beta = stof(argv[14]);
    }
    if(argc >= 14){
        visco0 = stof(argv[12]);
        Np_max = stof(argv[13]);
    }
    if(argc > 11 ) {
        phi_grid = stoi(argv[1]);
        theta_grid = stoi(argv[2]);
        diffusion_coeff = stof(argv[3]);
        shear_rate = stof(argv[4]);
        aspect_ratio = stof(argv[5]);
        dt = stof(argv[6]);
        max_time = stof(argv[7]);
        s1 = stoi(argv[8]);
        s2 = stoi(argv[9]);
        classes = stoi(argv[10]);
        what_todo = stoi(argv[11]);
    }
    else{ 
        cout << "INVALID INPUT!! \n";
        return 1;
    }
    cout << "#executing " << argv[0] << " " << phi_grid << " " << theta_grid << " " << diffusion_coeff
    << " " << shear_rate << " " << aspect_ratio << " " << dt << " " 
    << max_time << " " << s1 << " " << s2 << " " << classes << " " 
    << what_todo << " " << visco0 << " " << Np_max << " " << beta << " " << shear_rate_max <<  endl;
    if(what_todo){
        cout << "#selfcheck triggered" << endl; 
    }        
    Combo_class combined;
    Combo_class *user_data;
    combined.initialize(aspect_ratio, shear_rate, diffusion_coeff, classes, phi_grid, theta_grid, s1, s2, what_todo,
                        visco0, Np_max, beta, shear_rate_max);
    user_data=&combined;
    int system_size=(phi_grid+1)*theta_grid*classes;

    cout << "#SYSTEM SIZE " << system_size << endl;
    cout <<"#shearrate " << shear_rate << endl;
    cout <<"#s1: " << s1 << " s2: " << s2 << endl;
    
    ///initialize sundials stuff
    void * cvode_mem = nullptr;
	N_Vector y0 = nullptr;
    y0 = N_VNew_Serial(system_size);

    //Default orientation distribution 
    //Uniform distrtibution -->
    for(int i = 0; i != system_size; ++i){
		combined.set_uniform_state(NV_DATA_S(y0));
    }
    //<-- Uniform distribution
    //If what_todo==1, load the previous state
    if(what_todo == 1){ 
        combined.load_state(NV_DATA_S(y0)); 
    }
    cvode_mem = CVodeCreate(CV_ADAMS, CV_FUNCTIONAL);
    CVodeSetUserData(cvode_mem, user_data);
    /// Initialize the integrator memory
    CVodeInit(cvode_mem, wfunc, t0, y0);
    /// Set the error weights
    CVodeSStolerances(cvode_mem, 1.0e-8, 1.0e-8);  
    CVodeSetMaxNumSteps(cvode_mem,100000);
    ///LINEAR SOLVER MODULE
    ///specify dense solver
    CVSpgmr(cvode_mem, 0, 0);
    /// We need to run the solver for very tiny timestep to initialize it properly:
    CVode(cvode_mem, (1e-14), y0, &time, CV_ONE_STEP);

    int i = 0;
    double print_dt = 0.01, next_print_at = 0;
    double visco_0 = 0.001;
    while(time < max_time){
        i++;
        int flag;
        flag = CVode(cvode_mem, (time + dt), y0, &time, CV_NORMAL);
        if(time > next_print_at){
            ///TODO implement some printing function
            double visc_raw_til=combined.get_visc_raw(s1,s2); //luetaan visc_raw muuttujaan
            double visc_tot=visco0 + 2*visco0*visc_raw_til;  //Käytä: visc_raw Np_max visco0
            cout << "#AT time " << time << endl;
            cout << "#VISCORAW " << time << " " << combined.get_visc_raw(s1,s2) << endl;
            cout << "#VISCOTOT " << time << " " << visc_tot << endl;
            combined.save_aggr_distribution(time);
            next_print_at+=print_dt;
        }
    }
    combined.save_state();
    N_VDestroy_Serial(y0);
    CVodeFree(&cvode_mem);
    
    return 0;
}

//The function needed by the sundials solver "return dt"
//(get ydot at given t, y, and *userdata)
int wfunc (realtype t, N_Vector y, N_Vector ydot, void* userdata)
{
    auto * combined_ptr = (Combo_class*) userdata;
    combined_ptr->get_dt(NV_DATA_S(ydot), NV_DATA_S(y));
    return 0;
}