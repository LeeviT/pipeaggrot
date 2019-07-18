#include <iostream>
#include <cstdio>
#include <cmath>
#include <fstream>
#include <array>
#include <ctime>
#include <cstdlib>
#include <vector>

#include "point.hpp"
#include "io.hpp"
#include "velo.hpp"
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
int wfunc(realtype t, N_Vector y, N_Vector ydot, void * userdata) {
    auto * combined_ptr = (Combo_class *) userdata;
    combined_ptr->get_dt(NV_DATA_S(ydot), NV_DATA_S(y));
    return 0;
}

// Initialize viscosities containing vector which points radially from the midpoint outwards
// Viscosity is homogenous in the beginning, coming from the input file
vector<double> init_visc_vector(int ny, double input_visc) {
    vector<double> visc;

    for (int i = 0; i <= ny; i++) {
        visc.push_back(input_visc);
    }
    return visc;
}

int main() {
    // Dummy variable for point class objects
    double r;
    // Initialize variables and class instances
    input params{};
    vector<Combo_class> combined;
    vector<Combo_class *> user_data;
    vector<double> visc_vector;
    vector<point> radius;
    double time = 0, shear_rate = 0.0, visc_raw_til;
    int system_size;
    // Initialize SUNDIALS stuff
    vector<void *> cvode_mem;
    vector<N_Vector> y0;
    vector<ofstream> visctotfile;
    // Startup velocity stuff
    int n_bessel_zeros = 200;

    // Read input file and assign input values to input::params struct
    params.read_input_file();

    // Print input values to the screen
    printf("dp=%f, l=%f, visc=%f, R=%f, dt=%f, diff_coeff=%f\n", params.getdp(), params.getl(), params.getvisc0(),
            params.getR(), params.getdt(), params.getdiff_coeff());
    printf("aspect_ratio=%f, beta=%f, shear_rate_max=%f, Np_max=%f\n", params.getaspect_ratio(), params.getbeta(),
            params.getshear_rate_max(), params.getNp_max());
    printf("theta_grid=%i, phi_grid=%i, classes=%i, what_todo=%i\n", params.gettheta_grid(), params.getphi_grid(),
            params.getclasses(), params.getwhat_todo());
    printf("s1=%i, s2=%i, ny=%i, t0=%i, t_max=%i\n\n", params.gets1(), params.gets2(), params.getny(),
            params.gett0(), params.gett_max());

    // Allocate memory and initialize uniform viscosity vector for the first timestep
    visc_vector = init_visc_vector(params.getny() + 1, params.getvisc0());

    // Startup velocity stuff
    vector<double> bessel_zeros(n_bessel_zeros), J_1_values(n_bessel_zeros);
    vector<vector<double>> J_0_values(params.getny() + 1);
    for (int i = 0; i <= params.getny(); i++) { J_0_values[i].resize(n_bessel_zeros); }
    bessel_zeros = calc_bessel_zeros(n_bessel_zeros);
    J_0_values = J_0(params.getny(), bessel_zeros);
    J_1_values = J_1(bessel_zeros);

    // Initialize radius and combined vectors, ny+1 points along radius of the pipe
    // This must be done using *_back() function since vectors are accessed index-wise later on
    for (int i = 0; i <= params.getny(); i++) {
        radius.emplace_back();
        combined.emplace_back();
        visctotfile.emplace_back();
    }

    // Create ny AO-populations
    for (int i = 0; i <= params.getny(); i++) {
        // Population system size in spherical coordinates
        system_size = (params.getphi_grid() + 1) * params.gettheta_grid() * params.getclasses();
        // Initialize AO model class instance (single population) using input values
        combined.at(i).initialize(shear_rate, params);
        user_data.push_back(&combined.at(i));
        y0.push_back(N_VNew_Serial(system_size));
        // Set uniform distribution
        for (int j = 0; j < system_size; j++) { combined.at(i).set_uniform_state(NV_DATA_S(y0.at(i))); }
        cvode_mem.push_back(CVodeCreate(CV_ADAMS, CV_FUNCTIONAL));
        CVodeSetUserData(cvode_mem.at(i), user_data.at(i));
        // Initialize the integrator memory
        CVodeInit(cvode_mem.at(i), wfunc, params.gett0(), y0.at(i));
        // Set the error weights
        CVodeSStolerances(cvode_mem.at(i), 1.0e-8, 1.0e-8);
        CVodeSetMaxNumSteps(cvode_mem.at(i), 100000);
        // Specify dense solver
        CVSpgmr(cvode_mem.at(i), 0, 0);
        // We need to run the solver for very tiny timestep to initialize it properly:
        CVode(cvode_mem.at(i), 1e-14, y0.at(i), &time, CV_ONE_STEP);
    }

    // The main timestep loop
    for (int timestep = 0; timestep <= params.gett_max(); timestep++) {
        printf("t_step=%i of %i\n", timestep, params.gett_max());
        cout << "time: " << time << endl;

        // Loop goes through the other points in y-direction, starting from the middle and sets the following
        // values for each discretization point
        for (int i = 0; i <= params.getny(); i++) {
            // If in y-midpoint of the pipe, r coordinate is zero, naturally
            r = (i == 0) ? 0.0 : ((double) i / params.getny()) * params.getR();
            radius[i].setr(r);
            radius[i].setvisc(visc_vector[i]);
            radius[i].setx(1.0);
            radius[i].setvx(vx_pipe(radius[i], params, bessel_zeros, J_1_values, J_0_values, i, timestep));

            // Update shear rate for the AO model
            combined[i].update_shear_rate(radius[i].getvx());

            // Some AO model solver stuff
            if (i == 0) { CVode(cvode_mem[i], time + params.getdt(), y0[i], &time, CV_NORMAL); }
            else { CVode(cvode_mem[i], time, y0[i], &time, CV_NORMAL); }
            visc_raw_til = combined[i].get_visc_raw(params.gets1(), params.gets2(), timestep, i);

            // Compute total viscosity of the fluid
            visc_vector[i] = params.getvisc0() + 2 * params.getvisc0() * visc_raw_til;

            // Some printing stuff
            visctotfile[i].open("../pysrc/visctot" + to_string(i) + ".dat", ios::app);
            visctotfile[i] << time << "\t" << visc_vector[i] << endl;
            visctotfile[i].close();
            cout << "r: " << radius[i].getr() << ", total visc: " << visc_vector[i] << ", v: "
                 << radius[i].getvx() << endl;
            combined[i].save_aggr_distribution(time);
        }
        cout << " " << endl;
        // Writes r and vx values to file
        write_r_vx(params.getny(), timestep, radius);
    }

    for (int i = 0; i <= params.getny(); i++) {
        combined[i].save_state();
        N_VDestroy_Serial(y0[i]);
        CVodeFree(&cvode_mem[i]);
    }

    return 0;
}
