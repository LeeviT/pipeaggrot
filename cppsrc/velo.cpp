#include <utility>

#include <utility>

#include <utility>

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

#include <boost/math/special_functions/bessel.hpp>
using namespace std;


// Calcute m zeros of the first kind Bessel functions of order n=1
vector<double> calc_bessel_zeros(int n) {
    vector<double> zeros(n);

    // Note that cyl_bessel_j_zero() starts iterator indexing from 1, i.e. zeros[0] would cause segfault
    boost::math::cyl_bessel_j_zero<double>(0.0, 1, n, zeros.begin());
    return zeros;
}

// Calculate the Bessel function of the first kind of order zero values at points λ_n*(r/R)
vector<vector<double>> J_0(int ny, vector<double> bessel_zeros) {
    vector<vector<double>> values(ny+1);

    // First initialize J_0(λ) vectors for each r/R value
    for (int i = 0; i <= ny; i++) {
        values[i].resize(bessel_zeros.size());
    }
    // Now assign J_0(λ_n*(r/R)) values to J_0(λ) vectors
    for (int i = 0; i <= ny; i++) {
        for (int j = 0; j < bessel_zeros.size(); j++) {
            values[i][j] = boost::math::cyl_bessel_j<double>(0, bessel_zeros[j]*(((double) i)/((double) ny)));
        }
    }
    return values;
}

// Calculate the Bessel function of the first kind of order one values at points λ_n
vector<double> J_1(vector<double> bessel_zeros) {
    vector<double> values(bessel_zeros.size());

    // Calculate J_1(λ_n) values
    for (int i = 0; i < bessel_zeros.size(); i++) {
        values[i] = boost::math::cyl_bessel_j<double>(1, bessel_zeros[i]);
    }
    return values;
}

// Calculate value of Fourier-Bessel series for each r value
double fourier_bessel(vector<double> bessel_zeros, vector<double> J_1_values, vector<vector<double>> J_0_values,
                      int r_i, int t_step, input params, double visc) {
    double sum = 0.0, dens = 1.0, partsum = 0.0;

    for (int i = 0; i < bessel_zeros.size(); i++) {
        partsum += J_0_values[r_i][i] / (J_1_values[i]*pow(bessel_zeros[i], 3));
        sum += (1.0/pow(bessel_zeros[i], 3)) * (J_0_values[r_i][i] / J_1_values[i])
               * exp(-pow(bessel_zeros[i], 2)*((visc/dens) * (params.getdt()*(double) t_step))/pow(params.getR(), 2));
    }
    sum = 2*sum*params.getdp()*pow(params.getR(), 2)*params.getR()/(params.getl()*visc);

    return sum;
}

// Function to calculate fluid x-velocity in a single point in the pipe
double vx_pipe(point single_point, input params, vector<double> bessel_zeros, vector<double> J_1,
               vector<vector<double>> J_0, int r_i, int t_step) {

    double fb_sum = fourier_bessel(move(bessel_zeros), move(J_1), move(J_0), r_i, t_step,
                                   params, single_point.getvisc());
    return ((params.getdp()/(4* single_point.getvisc()*params.getl()))*(pow(params.getR(), 2) -
            pow(single_point.getr(), 2))) - fb_sum;
}
