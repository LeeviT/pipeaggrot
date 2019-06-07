#include <iostream>
#include <cstdio>
#include <cmath>
#include <fstream>
#include <array>
#include <ctime>
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

double calc_ao(int t_step, double visct) {
    // beta = v_d/v_c
    double t0 = 0, dt = 1e-1, time = 0, max_time = 1, print_dt = 0.01;
    double diffusion_coeff = 10, aspect_ratio = 20, shear_rate = 100, beta = 1.0, shear_rate_max = 100.0;
    double Np_max = 1.0, visc_tot, visc_raw_til, visc0 = 1.0, visc;
    int theta_grid = 20, phi_grid = 20, classes = 1, what_todo = 0, s1 = 0, s2 = 1;
    int system_size = (phi_grid + 1)*theta_grid*classes;

    if (t_step == 0) {
        visc = visc0;
    } else {
        visc = visct;
    }

    Combo_class combined;
    Combo_class *user_data;
    combined.initialize(aspect_ratio, shear_rate, diffusion_coeff, classes, phi_grid, theta_grid, s1, s2,
                        what_todo, visc, Np_max, beta, shear_rate_max);
    user_data = &combined;

    // initialize sundials stuff
    void * cvode_mem = nullptr;
    N_Vector y0 = nullptr;
    y0 = N_VNew_Serial(system_size);
    time += t_step*0.1;

    // Default orientation distribution
    // Set uniform distribution
    for (int i = 0; i != system_size; ++i) {
        combined.set_uniform_state(NV_DATA_S(y0));
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

    visc_raw_til = combined.get_visc_raw(s1, s2); //luetaan visc_raw muuttujaan
    visc_tot = visc + 2*visc*visc_raw_til;  //Käytä: visc_raw Np_max visc0
    cout << "#AT time " << time << endl;
    cout << "#VISCOTOT " << time << " " << visc_tot << endl;
    combined.save_aggr_distribution(time);

    combined.save_state();
    N_VDestroy_Serial(y0);
    CVodeFree(&cvode_mem);

    return visc_tot;
}

// Class to describe properties of a single discretized point in space (in cylindrical coord's)
// Contains spatial coordinates (x, r), viscosity (visc), velocity in x-direction (vx)
class point {
private:
    double x, r, visc, vx;

public:
    // Function for printing values of a point, mostly for debugging
    void print_values() const {
        printf("x=%f, r=%f, visc=%f, vx=%f\n", x, r, visc, vx);
    };
    // Setter functions for each variable
    void setx(double x_set) { x = x_set; };
    void setr(double r_set) { r = r_set; };
    void setvisc(double visc_set) { visc = visc_set; };
    void setvx(double vx_set) { vx = vx_set; }
    // Getter functions for each variable
    double getx() { return x; };
    double getr() { return r; };
    double getvisc() { return visc; };
    double getvx() { return vx; };
};


// Function to calculate fluid x-velocity in a single point in the pipe
double vx_pipe(point single_point, double dp, double l, double radius) {
    return (dp/(4*single_point.getvisc()*l))*(pow(radius, 2) - pow(single_point.getr(), 2));
}

// Writes r values and x-velocities along the radius to text file
void write_r_vx(int ny, int t_step, point radius[]) {
    ofstream file;
    // string filename = "../pysrc/r_vx" + to_string(t_step) + ".dat";
    file.open("../pysrc/r_vx" + to_string(t_step) + ".dat");
    for (int i = 0; i <= ny; i++) {
        file << radius[i].getr() << "\t" << radius[i].getvx() << "\n";
    }
    file.close();
}

// Structure for handling input and returning it
struct input {
    double dp, l, visc, R;
    int ny, t_max;
};

// Reads values from the input file to input structure and returns it
input read_input_file() {
    string a, b;
    double c;
    input values{};
    ifstream file;
    file.open("input.txt");

    while (file >> a >> b >> c) {
        if (a == "pressure_diff") { values.dp = c; }
        if (a == "pipe_length") { values.l = c; }
        if (a == "visc") { values.visc = c; }
        if (a == "radius") { values.R = c; }
        if (a == "ny") { values.ny = (int) c; }
        if (a == "t_max") { values.t_max = (int) c; }
    }
    file.close();
    return values;
}

// Initialize viscosities containing vector which points radially from the midpoint outwards
// Viscosity is homogenous, coming from the input file
double *init_visc_vector(int ny, double input_visc) {
    auto *visc = new double[ny + 1];

    for (int i = 0; i <= ny; i++) {
        visc[i] = input_visc;
    }
    return visc;
}

// Sets modified viscosity values to viscosity vector each timestep
// This will be the interface function to get viscosities from AO model
double *set_visc_vector(int t_step, int ny, double *visc_vector) {

    for (int i = 0; i <= ny; i++) {
        visc_vector[i] = calc_ao(t_step, visc_vector[i]);
    }
    return visc_vector;
}

int main() {
    // Variables for point class objects
    double r;
    // Initial values and constants for the system from the input file
    // dp=pressure difference, l=length of the pipe, R=radius of the pipe
    // nx=number of discretization points in x-direction, ny=same but in y-direction, t_max=max number of timesteps
    input input_params{};
    double dp, l, visc, R;
    double *visc_vector;
    int ny, t_max;

    // Assigns input values to local variables using input structure
    // Some values won't change during program runs so maybe could use values straight from struct, more error prone?
    input_params = read_input_file();
    dp = input_params.dp;
    l = input_params.l;
    visc = input_params.visc;
    R = input_params.R;
    ny = input_params.ny;
    t_max = input_params.t_max;
    printf("dp=%f, l=%f, visc=%f, R=%f, ny=%i, t_max=%i\n", dp, l, visc, R, ny, t_max);

    // Allocate memory and initialize viscosity vector for the first timestep
    visc_vector = init_visc_vector(ny + 1, visc);

    // Create point object, ny points along radius of the pipe
    point radius[ny + 1];

    for (int t_step = 0; t_step <= t_max; t_step++) {
        printf("t_step=%i of %i\n", t_step, t_max);
        // Set new viscosity values if not in the very first timestep
        if (t_step > 0) {
            set_visc_vector(t_step, ny, visc_vector);
        }

        // Sets values for the point in the y-midpoint
        radius[0].setr(0.0);
        radius[0].setvisc(visc_vector[0]);
        radius[0].setx(1.0);
        radius[0].setvx(vx_pipe(radius[0], dp, l, R));

        // Loop goes through the other points in y-direction, starting from the middle
        for (int i = 1; i <= ny; i++) {
            r = ((double) i / ny) * R;
            radius[i].setr(r);
            radius[i].setvisc(visc_vector[i]);
            radius[i].setx(1.0);
            radius[i].setvx(vx_pipe(radius[i], dp, l, R));
        }

        // Prints values along the radius
        for (int i = 0; i <= ny; i++) {
            // radius[i].print_values();
        }
        // Writes r and values to file
        write_r_vx(ny, t_step, radius);
    }

    // delete [] visc_vector;

    return 0;
}
