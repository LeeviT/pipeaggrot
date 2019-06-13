#include <iostream>
#include <cstdio>
#include <cmath>
#include <fstream>
#include <array>
#include <ctime>
#include <cstdlib>
#include <vector>

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
    void setx(double x) { point::x = x; };
    void setr(double r) { point::r = r; };
    void setvisc(double visc) { point::visc = visc; };
    void setvx(double vx) { point::vx = vx; }
    // Getter functions for each variable
    double getx() const { return x; };
    double getr() const { return r; };
    double getvisc() const { return visc; };
    double getvx() const { return vx; };
};

// Structure for handling input and returning it
struct input {
    // beta = v_d/v_c
    // dp=pressure difference, l=length of the pipe, R=radius of the pipe
    // nx=number of discretization points in x-direction, ny=same but in y-direction, t_max=max number of timesteps
    double dp, l, visc0, R, dt, diff_coeff, aspect_ratio, beta, shear_rate_max, Np_max;
    int theta_grid, phi_grid, classes, what_todo, s1, s2, ny, t0, t_max;

    void setdp(double dp) { input::dp = dp; }
    void setl(double l) { input::l = l; }
    void setvisc0(double visc0) { input::visc0 = visc0; }
    void setR(double R) { input::R = R; }
    void setdt(double dt) { input::dt = dt; }
    void setdiff_coeff(double diff_coeff) { input::diff_coeff = diff_coeff; }
    void setaspect_ratio(double aspect_ratio) { input::aspect_ratio = aspect_ratio; }
    void setbeta(double beta) { input::beta = beta;}
    void setshear_rate_max(double shear_rate_max) { input::shear_rate_max = shear_rate_max; }
    void setNp_max(double Np_max) { input::Np_max = Np_max; }
    void settheta_grid(int theta_grid) { input::theta_grid = theta_grid; }
    void setphi_grid(int phi_grid) { input::phi_grid = phi_grid; }
    void setclasses(int classes) { input::classes = classes; }
    void setwhat_todo(int what_todo) { input::what_todo = what_todo;}
    void sets1(int s1) { input::s1 = s1; }
    void sets2(int s2) { input::s2 = s2; }
    void setny(int ny) { input::ny = ny; }
    void sett0(int t0) { input::t0 = t0; }
    void sett_max(int t_max) { input::t_max = t_max; }

    double getdp() const { return dp; }
    double getl() const { return l; }
    double getvisc0() const { return visc0; }
    double getR() const { return R; }
    double getdt() const { return dt; }
    double getdiff_coeff() const { return diff_coeff; }
    double getaspect_ratio() const { return aspect_ratio; }
    double getbeta() const { return beta; }
    double getshear_rate_max() const { return shear_rate_max; }
    double getNp_max() const { return Np_max; }
    int gettheta_grid() const { return theta_grid; }
    int getphi_grid() const { return phi_grid; }
    int getclasses() const { return classes; }
    int getwhat_todo() const { return what_todo; }
    int gets1() const { return s1; }
    int gets2() const { return s2; }
    int getny() const { return ny; }
    int gett0() const { return t0; }
    int gett_max() const { return t_max; }
};

// Function to calculate fluid x-velocity in a single point in the pipe
double vx_pipe(point single_point, double dp, double l, double radius) {
    return (dp/(4* single_point.getvisc()*l))*(pow(radius, 2) - pow(single_point.getr(), 2));
}

// Writes r values and x-velocities along the radius to text file
void write_r_vx(int ny, int t_step, point radius[]) {
    ofstream file;
    file.open("../pysrc/r_vx" + to_string(t_step) + ".dat");
    for (int i = 0; i <= ny; i++) {
        file << radius[i].getr() << "\t" << radius[i].getvx() << "\n";
    }
    file.close();
}

// Reads values from the input file to input structure and returns it
input read_input_file() {
    string a, b;
    double c;
    input values{};
    ifstream file;
    file.open("input.txt");

    while (file >> a >> b >> c) {
        if (a == "pressure_diff") { values.setdp(c); }
        if (a == "pipe_length") { values.setl(c); }
        if (a == "visc") { values.setvisc0(c); }
        if (a == "radius") { values.setR(c); }
        if (a == "ny") { values.setny((int) c); }
        if (a == "t_max") { values.sett_max((int) c); }
    }
    file.close();
    return values;
}

// Initialize viscosities containing vector which points radially from the midpoint outwards
// Viscosity is homogenous, coming from the input file
vector<double> init_visc_vector(int ny, double input_visc) {
    vector<double> visc;

    for (int i = 0; i <= ny; i++) {
        visc.push_back(input_visc);
    }
    return visc;
}

// Sets modified viscosity values to viscosity vector each timestep
// This will be the interface function to get viscosities from AO model
vector<double> set_visc_vector(int ny, double visc, vector<double> visc_vector) {

    for (int i = 0; i <= ny; i++) {
        visc_vector[i] = visc;
    }
    return visc_vector;
}

int main() {
    // Variables for point class objects
    double r;
    // Initial values and constants for the system from the input file
    ofstream visctotfile;
    Combo_class combined;
    Combo_class *user_data;
    input params{};
    vector<double> visc_vector;
    double time = 0, shear_rate = 100.0, visc_raw_til, visc;
    int system_size;

    // Assigns input values to local variables using input structure
    // Some values won't change during program runs so maybe could use values straight from struct, more error prone?
    params = read_input_file();
    printf("dp=%f, l=%f, visc=%f, R=%f, ny=%i, t_max=%i\n", params.getdp(), params.getl(), params.getvisc0(),
            params.getR(), params.getny(), params.gett_max());

    // Allocate memory and initialize viscosity vector for the first timestep
    visc_vector = init_visc_vector(params.getny() + 1, params.getvisc0());

    // Create point object, ny points along radius of the pipe
    point radius[params.getny() + 1];

    visc = params.getvisc0();
    system_size = (params.getphi_grid() + 1)*params.gettheta_grid()*params.getclasses();
    combined.initialize(params.getaspect_ratio(), shear_rate, params.getdiff_coeff(), params.getclasses(),
                        params.getphi_grid(), params.gettheta_grid(), params.gets1(), params.gets2(),
                        params.getwhat_todo(), params.getvisc0(), params.getNp_max(), params.getbeta(),
                        params.getshear_rate_max());
    user_data = &combined;
    // initialize sundials stuff
    void * cvode_mem = nullptr;
    N_Vector y0 = nullptr;
    y0 = N_VNew_Serial(system_size);
    // Set uniform distribution
    for (int i = 0; i != system_size; ++i) {
        combined.set_uniform_state(NV_DATA_S(y0));
    }
    cvode_mem = CVodeCreate(CV_ADAMS, CV_FUNCTIONAL);
    CVodeSetUserData(cvode_mem, user_data);
    /// Initialize the integrator memory
    CVodeInit(cvode_mem, wfunc, params.gett0(), y0);
    /// Set the error weights
    CVodeSStolerances(cvode_mem, 1.0e-8, 1.0e-8);
    CVodeSetMaxNumSteps(cvode_mem, 100000);
    ///specify dense solver
    CVSpgmr(cvode_mem, 0, 0);
    /// We need to run the solver for very tiny timestep to initialize it properly:
    CVode(cvode_mem, 1e-14, y0, &time, CV_ONE_STEP);

    // Open file to write visctot values
    visctotfile.open("../pysrc/visctot.dat");

    // The main timestep loop
    for (int t_step = 0; t_step <= params.gett_max(); t_step++) {
        printf("t_step=%i of %i\n", t_step, params.gett_max());

        // Set new viscosity values if not in the very first timestep
        if (t_step > 0) {
            set_visc_vector(params.getny(), visc, visc_vector);
        }

        // Sets values for the point in the y-midpoint
        radius[0].setr(0.0);
        radius[0].setvisc(visc_vector[0]);
        radius[0].setx(1.0);
        radius[0].setvx(vx_pipe(radius[0], params.getdp(), params.getl(), params.getR()));

        // Loop goes through the other points in y-direction, starting from the middle
        for (int i = 1; i <= params.getny(); i++) {
            r = ((double) i / params.getny()) * params.getR();
            radius[i].setr(r);
            radius[i].setvisc(visc_vector[i]);
            radius[i].setx(1.0);
            radius[i].setvx(vx_pipe(radius[i], params.getdp(), params.getl(), params.getR()));
        }

        CVode(cvode_mem, (time + params.getdt()), y0, &time, CV_NORMAL);
        visc_raw_til = combined.get_visc_raw(params.gets1(), params.gets2()); //luetaan visc_raw muuttujaan
        visc = params.getvisc0() + 2*params.getvisc0()*visc_raw_til;  //Käytä: visc_raw Np_max visc0
        visctotfile << time << "\t" <<visc << endl;
        cout << "#AT time " << time << endl;
        cout << "#VISCOTOT " << time << " " << visc << endl;
        combined.save_aggr_distribution(time);

        // Writes r and vx values to file
        write_r_vx(params.getny(), t_step, radius);
    }
    visctotfile.close();

    combined.save_state();
    N_VDestroy_Serial(y0);
    CVodeFree(&cvode_mem);

    return 0;
}