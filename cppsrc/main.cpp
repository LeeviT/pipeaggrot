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
    // Constructor
    point() {
        point::x = 0.0;
        point::r = 0.0;
        point::visc = 0.0;
        point::vx = 0.0;
    }
    // Function for printing values of a point, mostly for debugging
    void print_values() const { printf("x=%f, r=%f, visc=%f, vx=%f\n", x, r, visc, vx); };
    // Setter functions for each variable
    void setx(double x) { point::x = x; };
    void setr(double r) { point::r = r; };
    void setvisc(double visc) { point::visc = visc; };
    void setvx(double vx) { point::vx = vx; }
    // Getter functions for each variable
    double getx() { return x; };
    double getr() { return r; };
    double getvisc() { return visc; };
    double getvx() { return vx; };
};

// Structure for handling input and returning it
struct input {
    // beta = v_d/v_c
    // dp=pressure difference, l=length of the pipe, R=radius of the pipe
    // nx=number of discretization points in x-direction, ny=same but in y-direction, t_max=max number of timesteps
    double dp, l, visc0, R, dt, diff_coeff, aspect_ratio, beta, shear_rate_max, Np_max;
    int theta_grid, phi_grid, classes, what_todo, s1, s2, ny, t0, t_max;
    // Setter functions for each variable
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
    // Getter functions for each variable
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
void write_r_vx(int ny, int t_step, vector<point> radius) {
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

    // Goes through the input file line by line and breaks every line into "variable name", "=" and actual
    // input value, then passes input values into input::values struct
    while (file >> a >> b >> c) {
        if (a == "pressure_diff") { values.setdp(c); }
        if (a == "pipe_length") { values.setl(c); }
        if (a == "visc") { values.setvisc0(c); }
        if (a == "radius") { values.setR(c); }
        if (a == "dt") { values.setdt(c); }
        if (a == "diffusion_coeff") { values.setdiff_coeff(c); }
        if (a == "aspect_ratio") { values.setaspect_ratio(c); }
        if (a == "beta") { values.setbeta(c); }
        if (a == "shear_rate_max") { values.setshear_rate_max(c); }
        if (a == "Np_max") { values.setNp_max(c); }
        if (a == "ny") { values.setny((int) c); }
        if (a == "t_max") { values.sett_max((int) c); }
        if (a == "theta_grid") { values.settheta_grid((int) c); }
        if (a == "phi_grid") { values.setphi_grid((int) c); }
        if (a == "classes") { values.setclasses((int) c); }
        if (a == "what_todo") { values.setwhat_todo((int) c); }
        if (a == "s1") { values.sets1((int) c); }
        if (a == "s2") { values.sets2((int) c); }
        if (a == "t0") { values.sett0((int) c); }
    }
    file.close();
    return values;
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

// Sets modified viscosity values to viscosity vector each timestep
// This will be the interface function to get viscosities from AO model
vector<double> set_visc_vector(int ny, double visc) {
    vector<double> visc_vector;

    for (int i = 0; i <= ny; i++) {
        visc_vector.push_back(visc);
    }
    return visc_vector;
}

int main() {
    // Dummy variable for point class objects
    double r;
    // Initialize variables and class instances
    ofstream visctotfile;
    input params{};
    vector<Combo_class> combined;
    vector<Combo_class *> user_data;
    vector<double> visc_vector;
    vector<point> radius;
    double time = 0, shear_rate = 100.0, visc_raw_til, visc;
    int system_size;
    // Initialize sundials stuff
    vector<void *> cvode_mem;
    vector<N_Vector> y0;

    // Read input file and assign input values to input::params struct
    params = read_input_file();
    // Print input values to the screen
    printf("dp=%f, l=%f, visc=%f, R=%f, dt=%f, diff_coeff=%f\n", params.getdp(), params.getl(), params.getvisc0(),
            params.getR(), params.getdt(), params.getdiff_coeff());
    printf("aspect_ratio=%f, beta=%f, shear_rate_max=%f, Np_max=%f\n", params.getaspect_ratio(), params.getbeta(),
            params.getshear_rate_max(), params.getNp_max());
    printf("theta_grid=%i, phi_grid=%i, classes=%i, what_todo=%i\n", params.gettheta_grid(), params.getphi_grid(),
            params.getclasses(), params.getwhat_todo());
    printf("s1=%i, s2=%i, ny=%i, t0=%i, t_max=%i\n", params.gets1(), params.gets2(), params.getny(), params.gett0(),
            params.gett_max());

    // Allocate memory and initialize viscosity vector for the first timestep
    visc_vector = init_visc_vector(params.getny() + 1, params.getvisc0());

    // Initialize radius and combined vectors, ny+1 points along radius of the pipe
    // This must be done using *_back() function since radius vector is accessed index-wise later on
    for (int i = 0; i <= params.getny(); i++) {
        radius.emplace_back();
        combined.emplace_back();
    }

    // Initialize viscosity value based on the input viscosity
    visc = params.getvisc0();

    // Create ny population
    for (int i = 0; i <= params.getny(); i++) {
        // Population system size in spherical coordinates
        system_size = (params.getphi_grid() + 1) * params.gettheta_grid() * params.getclasses();
        // Initialize AO model class instance (single population) using input values
        combined.at(i).initialize(params.getaspect_ratio(), shear_rate, params.getdiff_coeff(), params.getclasses(),
                            params.getphi_grid(), params.gettheta_grid(), params.gets1(), params.gets2(),
                            params.getwhat_todo(), params.getvisc0(), params.getNp_max(), params.getbeta(),
                            params.getshear_rate_max());
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

    // Open file to write visctot values
    visctotfile.open("../pysrc/visctot.dat");

    // The main timestep loop
    for (int t_step = 0; t_step <= params.gett_max(); t_step++) {
        printf("t_step=%i of %i\n", t_step, params.gett_max());

        // Set new viscosity values if not in the very first timestep
        if (t_step > 0) { visc_vector = set_visc_vector(params.getny(), visc); }

        // Loop goes through the other points in y-direction, starting from the middle and sets the following
        // values for each discretization point
        for (int i = 0; i <= params.getny(); i++) {
            // If in y-midpoint of the pipe, r coordinate is zero, naturally
            if (i == 0) {
                r = 0.0;
            } else {
                r = ((double) i / params.getny()) * params.getR();
            }
            radius[i].setr(r);
            radius[i].setvisc(visc_vector[i]);
            radius[i].setx(1.0);
            radius[i].setvx(vx_pipe(radius[i], params.getdp(), params.getl(), params.getR()));

            // Update shear rate for the AO model
            combined.at(i).update_shear_rate(radius[i].getvx());

            // Some AO model solver stuff
            CVode(cvode_mem.at(i), (time + params.getdt()), y0.at(i), &time, CV_NORMAL);
            visc_raw_til = combined.at(i).get_visc_raw(params.gets1(), params.gets2()); //luetaan visc_raw muuttujaan

            // Compute total viscosity of the fluid
            visc = params.getvisc0() + 2 * params.getvisc0() * visc_raw_til;  //Käytä: visc_raw Np_max visc0

            // Some printing stuff
            visctotfile << time << "\t" << visc << endl;
            cout << "#AT time " << time << endl;
            cout << "#VISCOTOT " << time << " " << visc << endl;
            combined.at(i).save_aggr_distribution(time);
        }

        // Writes r and vx values to file
        write_r_vx(params.getny(), t_step, radius);
    }
    visctotfile.close();

    for (int i = 0; i <= params.getny(); i++) {
        combined.at(i).save_state();
        N_VDestroy_Serial(y0.at(i));
        CVodeFree(&cvode_mem.at(i));
    }

    return 0;
}
