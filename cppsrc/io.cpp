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
void input::setbeta(double beta_) { beta = beta_;}
void input::setshear_rate_max(double shear_rate_max) { input::shear_rate_max = shear_rate_max; }
void input::setNp_max(double Np_max) { input::Np_max = Np_max; }
void input::settheta_grid(int theta_grid) { input::theta_grid = theta_grid; }
void input::setphi_grid(int phi_grid) { input::phi_grid = phi_grid; }
void input::setclasses(int classes) { input::classes = classes; }
void input::setwhat_todo(int what_todo) { input::what_todo = what_todo;}
void input::sets1(int s1) { input::s1 = s1; }
void input::sets2(int s2) { input::s2 = s2; }
void input::setny(int ny) { input::ny = ny; }
void input::sett0(int t0) { input::t0 = t0; }
void input::sett_max(int t_max) { input::t_max = t_max; }

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

// Reads values from the input file to input structure and returns it
void input::read_input_file() {
    string a, b;
    double c;
    ifstream file;
    file.open("input.txt");

    // Goes through the input file line by line and breaks every line into "variable name", "=" and actual
    // input value, then passes input values into input::values struct
    while (file >> a >> b >> c) {
        if (a == "pressure_diff") { setdp(c); }
        if (a == "pipe_length") { setl(c); }
        if (a == "visc") { setvisc0(c); }
        if (a == "radius") { setR(c); }
        if (a == "dt") { setdt(c); }
        if (a == "diffusion_coeff") { setdiff_coeff(c); }
        if (a == "aspect_ratio") { setaspect_ratio(c); }
        if (a == "beta") { setbeta(c); }
        if (a == "shear_rate_max") { setshear_rate_max(c); }
        if (a == "Np_max") { setNp_max(c); }
        if (a == "ny") { setny((int) c); }
        if (a == "t_max") { sett_max((int) c); }
        if (a == "theta_grid") { settheta_grid((int) c); }
        if (a == "phi_grid") { setphi_grid((int) c); }
        if (a == "classes") { setclasses((int) c); }
        if (a == "what_todo") { setwhat_todo((int) c); }
        if (a == "s1") { sets1((int) c); }
        if (a == "s2") { sets2((int) c); }
        if (a == "t0") { sett0((int) c); }
    }
    file.close();
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