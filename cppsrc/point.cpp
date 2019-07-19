#include <iostream>
#include <cstdio>
#include <cmath>
#include <fstream>
#include <array>
#include <ctime>
#include <cstdlib>
#include <vector>

#include "point.hpp"

using namespace std;

// Constructor
point::point() {
    x = 0.0;
    r = 0.0;
    visc = 0.0;
    vx = 0.0;
}

// Function for printing values of a point, mostly for debugging
void point::print_values() { printf("x=%f, r=%f, visc=%f, vx=%f\n", x, r, visc, vx); }

// Setter functions for each variable
void point::setx(double x_) { x = x_; }
void point::setr(double r_) { r = r_; }
void point::setvisc(double visc_) { visc = visc_; }
void point::setvx(double vx_) { vx = vx_; }

// Getter functions for each variable
double point::getx() { return x; }
double point::getr() { return r; }
double point::getvisc() { return visc; }
double point::getvx() { return vx; }
