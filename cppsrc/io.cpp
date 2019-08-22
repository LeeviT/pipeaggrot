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

using namespace std;

// Setter functions for each variable
void input::setdp(double dp_) { dp = dp_; }
void input::setl(double l_) { l = l_; }
void input::setvisc0(double visc0_) { visc0 = visc0_; }
void input::setR(double R_) { R = R_; }
void input::setdt(double dt_) { dt = dt_; }
void input::setdiff_coeff(double diff_coeff_) { diff_coeff = diff_coeff_; }
void input::setaspect_ratio(double aspect_ratio_) { aspect_ratio = aspect_ratio_; }
void input::setbeta(double beta_) { beta = beta_; }
void input::setshear_rate_max(double shear_rate_max_) { shear_rate_max = shear_rate_max_; }
void input::setNp_max(double Np_max_) { Np_max = Np_max_; }
void input::settheta_grid(int theta_grid_) { theta_grid = theta_grid_; }
void input::setphi_grid(int phi_grid_) { phi_grid = phi_grid_; }
void input::setclasses(int classes_) { classes = classes_; }
void input::setwhat_todo(int what_todo_) { what_todo = what_todo_; }
void input::sets1(int s1_) { s1 = s1_; }
void input::sets2(int s2_) { s2 = s2_; }
void input::setny(int ny_) { ny = ny_; }
void input::sett0(int t0_) { t0 = t0_; }
void input::sett_max(int t_max_) { t_max = t_max_; }

// Getter functions for each variable
double input::getdp() const { return dp; }
double input::getl() const { return l; }
double input::getvisc0() const { return visc0; }
double input::getR() const { return R; }
double input::getdt() const { return dt; }
double input::getdiff_coeff() const { return diff_coeff; }
double input::getaspect_ratio() const { return aspect_ratio; }
double input::getbeta() const { return beta; }
double input::getshear_rate_max() const { return shear_rate_max; }
double input::getNp_max() const { return Np_max; }
int input::gettheta_grid() const { return theta_grid; }
int input::getphi_grid() const { return phi_grid; }
int input::getclasses() const { return classes; }
int input::getwhat_todo() const { return what_todo; }
int input::gets1() const { return s1; }
int input::gets2() const { return s2; }
int input::getny() const { return ny; }
int input::gett0() const { return t0; }
int input::gett_max() const { return t_max; }

// Reads values from the input _file to input structure and returns it
void input::read_input_file() {
    string _a, _b;
    double _c;
    ifstream _file;

    _file.open("input.txt");

    // Goes through the input _file line by line and breaks every line into "variable name", "=" and actual
    // input value, then passes input values into input::values struct
    while (_file >> _a >> _b >> _c) {
        if (_a == "pressure_diff") { setdp(_c); }
        if (_a == "pipe_length") { setl(_c); }
        if (_a == "visc") { setvisc0(_c); }
        if (_a == "radius") { setR(_c); }
        if (_a == "dt") { setdt(_c); }
        if (_a == "diffusion_coeff") { setdiff_coeff(_c); }
        if (_a == "aspect_ratio") { setaspect_ratio(_c); }
        if (_a == "beta") { setbeta(_c); }
        if (_a == "shear_rate_max") { setshear_rate_max(_c); }
        if (_a == "Np_max") { setNp_max(_c); }
        if (_a == "ny") { setny((int) _c); }
        if (_a == "t_max") { sett_max((int) _c); }
        if (_a == "theta_grid") { settheta_grid((int) _c); }
        if (_a == "phi_grid") { setphi_grid((int) _c); }
        if (_a == "classes") { setclasses((int) _c); }
        if (_a == "what_todo") { setwhat_todo((int) _c); }
        if (_a == "s1") { sets1((int) _c); }
        if (_a == "s2") { sets2((int) _c); }
        if (_a == "t0") { sett0((int) _c); }
    }
    _file.close();
}

// Writes r values and x-velocities along the radius to text _file
void write_r_vx(int ny_, int t_step_, vector<point> radius_) {
    ofstream _file;

    _file.open("output/r_vx" + to_string(t_step_) + ".dat");

    for (int i = 0; i <= ny_; i++) {
        _file << radius_[i].getr() << "\t" << radius_[i].getvx() << "\n";
    }
    _file.close();
}