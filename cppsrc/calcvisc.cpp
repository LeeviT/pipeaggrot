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

// The function needed by the sundials solver "return dt"
// (get ydot at given t, y, and *userdata)
int wfunc (realtype t, N_Vector y, N_Vector ydot, void* userdata) {
    auto * combined_ptr = (Combo_class*) userdata;
    combined_ptr->get_dt(NV_DATA_S(ydot), NV_DATA_S(y));
    return 0;
}

double calc_ao() {
    // beta = v_d/v_c
    double t0 = 0, dt = 1e-1, time = 0, max_time = 1, print_dt = 0.01;
    double diffusion_coeff = 10, aspect_ratio = 20, shear_rate = 100, beta = 1.0, shear_rate_max = 100.0, Np_max = 1.0;
    double visc_tot, visc_raw_til, visco0 = 1.0;
    int theta_grid = 20, phi_grid = 20, classes = 1, what_todo = 0, s1 = 0, s2 = 1;
    int system_size = (phi_grid + 1)*theta_grid*classes;

    cout << "###########################################################\n";
    cout << "############ AGGREGATION ORIENTATION MODEL ################\n";
    cout << "###########################################################\n";

    cout << "#executing " << " " << phi_grid << " " << theta_grid << " " << diffusion_coeff << " " << shear_rate <<
    " " << aspect_ratio << " " << dt << " " << max_time << " " << s1 << " " << s2 << " " << classes <<
    " " << what_todo << " " << visco0 << " " << Np_max << " " << beta << " " << shear_rate_max <<  endl;

    Combo_class combined;
    Combo_class *user_data;
    combined.initialize(aspect_ratio, shear_rate, diffusion_coeff, classes, phi_grid, theta_grid, s1, s2,
                        what_todo, visco0, Np_max, beta, shear_rate_max);
    user_data = &combined;

    cout << "#SYSTEM SIZE " << system_size << endl;
    cout << "#shearrate " << shear_rate << endl;
    cout << "#s1: " << s1 << " s2: " << s2 << endl;
    
    // initialize sundials stuff
    void * cvode_mem = nullptr;
	N_Vector y0 = nullptr;
    y0 = N_VNew_Serial(system_size);

    // Default orientation distribution
    // Uniform distrtibution -->
    for (int i = 0; i != system_size; ++i) {
		combined.set_uniform_state(NV_DATA_S(y0));
    }
    // <-- Uniform distribution
    // If what_todo==1, load the previous state
    if (what_todo == 1) {
        combined.load_state(NV_DATA_S(y0));
    }
    cvode_mem = CVodeCreate(CV_ADAMS, CV_FUNCTIONAL);
    CVodeSetUserData(cvode_mem, user_data);
    /// Initialize the integrator memory
    CVodeInit(cvode_mem, wfunc, t0, y0);
    /// Set the error weights
    CVodeSStolerances(cvode_mem, 1.0e-8, 1.0e-8);  
    CVodeSetMaxNumSteps(cvode_mem, 100000);
    ///LINEAR SOLVER MODULE
    ///specify dense solver
    CVSpgmr(cvode_mem, 0, 0);
    /// We need to run the solver for very tiny timestep to initialize it properly:
    CVode(cvode_mem, 1e-14, y0, &time, CV_ONE_STEP);
    CVode(cvode_mem, (time + dt), y0, &time, CV_NORMAL);

    ///TODO implement some printing function
    visc_raw_til = combined.get_visc_raw(s1, s2); //luetaan visc_raw muuttujaan
    visc_tot = visco0 + 2*visco0*visc_raw_til;  //Käytä: visc_raw Np_max visc0
    cout << "#AT time " << time << endl;
    cout << "#VISCORAW " << time << " " << combined.get_visc_raw(s1, s2) << endl;
    cout << "#VISCOTOT " << time << " " << visc_tot << endl;
    combined.save_aggr_distribution(time);

    combined.save_state();
    N_VDestroy_Serial(y0);
    CVodeFree(&cvode_mem);

    return visc_tot;
}
