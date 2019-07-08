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
vector<double> calc_bessel_zeros(int n);

// Calculate the Bessel function of the first kind of order zero values at points λ_n*(r/R)
vector<vector<double>> J_0(int ny, vector<double> bessel_zeros);

// Calculate the Bessel function of the first kind of order one values at points λ_n
vector<double> J_1(vector<double> bessel_zeros);

// Calculate value of Fourier-Bessel series for each r value
double fourier_bessel(vector<double> bessel_zeros, vector<double> J_1_values, vector<vector<double>> J_0_values,
                      int r_i, int t_step, input params, double visc);

// Function to calculate fluid x-velocity in a single point in the pipe
double vx_pipe(point single_point, input params, vector<double> bessel_zeros, vector<double> J_1,
               vector<vector<double>> J_0, int r_i, int t_step);

#endif