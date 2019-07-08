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
    point::x = 0.0;
    point::r = 0.0;
    point::visc = 0.0;
    point::vx = 0.0;
}

// Function for printing values of a point, mostly for debugging
void point::print_values() { printf("x=%f, r=%f, visc=%f, vx=%f\n", x, r, visc, vx); }

// Setter functions for each variable
void point::setx(double x) { point::x = x; }
void point::setr(double r) { point::r = r; }
void point::setvisc(double visc) { point::visc = visc; }
void point::setvx(double vx) { point::vx = vx; }

// Getter functions for each variable
double point::getx() { return point::x; }
double point::getr() { return point::r; }
double point::getvisc() { return point::visc; }
double point::getvx() { return point::vx; }
