#ifndef CPPSRC_VELO_HPP
#define CPPSRC_VELO_HPP

#include <iostream>
#include <cstdio>
#include <cmath>
#include <fstream>
#include <array>
#include <ctime>
#include <cstdlib>
#include <vector>

using namespace std;

// cpp/hpp files containing functions to calculate startup velocity and stable state velocity of pipe flow
// Bessel function-related functions are used to calculate startup velocity

// Calcute m zeros of the first kind Bessel functions of order n=1
vector<double> calc_bessel_zeros(int n_);

// Calculate the Bessel function of the first kind of order zero values at points λ_n*(r/R)
vector<vector<double>> J_0(int ny_, vector<double> bessel_zeros_);

// Calculate the Bessel function of the first kind of order one values at points λ_n
vector<double> J_1(vector<double> bessel_zeros_);

// Calculate value of Fourier-Bessel series for each r value
double fourier_bessel(vector<double> bessel_zeros_, vector<double> J_1_values_, vector<vector<double>> J_0_values_,
                      int r_i_, int t_step_, input params_, double visc_);

// Function to calculate fluid x-velocity in a single point in the pipe
double vx_pipe(point single_point_, input params_, vector<double> bessel_zeros_, vector<double> J_1_,
               vector<vector<double>> J_0_, int r_i_, int t_step_);

#endif