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
vector<double> calc_bessel_zeros(int n_) {
    vector<double> _zeros(n_);

    // Note that cyl_bessel_j_zero() starts iterator indexing from 1, i.e. zeros[0] would cause segfault
    boost::math::cyl_bessel_j_zero<double>(0.0, 1, n_, _zeros.begin());
    return _zeros;
}

// Calculate the Bessel function of the first kind of order zero values at points λ_n*(r/R)
vector<vector<double>> J_0(int ny_, vector<double> bessel_zeros_) {
    vector<vector<double>> _values(ny_ + 1);

    // First initialize J_0(λ) vectors for each r/R value
    for (int i = 0; i <= ny_; i++) {
        _values[i].resize(bessel_zeros_.size());
    }
    // Now assign J_0(λ_n*(r/R)) values to J_0(λ) vectors
    for (int i = 0; i <= ny_; i++) {
        for (int j = 0; j < bessel_zeros_.size(); j++) {
            _values[i][j] = boost::math::cyl_bessel_j<double>(0, bessel_zeros_[j]*(((double) i)/((double) ny_)));
        }
    }
    return _values;
}

// Calculate the Bessel function of the first kind of order one values at points λ_n
vector<double> J_1(vector<double> bessel_zeros_) {
    vector<double> _values(bessel_zeros_.size());

    // Calculate J_1(λ_n) values
    for (int i = 0; i < bessel_zeros_.size(); i++) {
        _values[i] = boost::math::cyl_bessel_j<double>(1, bessel_zeros_[i]);
    }
    return _values;
}

// Calculate value of Fourier-Bessel series for each r value
double fourier_bessel(vector<double> bessel_zeros_, vector<double> J_1_values_, vector<vector<double>> J_0_values_,
                      int r_i_, int t_step_, input params_, double visc_) {
    double _sum = 0.0, _dens = 1.0, _partsum = 0.0;

    for (int i = 0; i < bessel_zeros_.size(); i++) {
        _partsum += J_0_values_[r_i_][i] / (J_1_values_[i] * pow(bessel_zeros_[i], 3));
        _sum += (1.0/pow(bessel_zeros_[i], 3)) * (J_0_values_[r_i_][i] / J_1_values_[i])
              * exp(-pow(bessel_zeros_[i], 2) *((visc_ / _dens) * (params_.getdt()*(double) t_step_))
              / pow(params_.getR(), 2));
    }
    _sum = 2 * _sum * params_.getdp() * pow(params_.getR(), 2) * params_.getR() / (params_.getl() * visc_);

    return _sum;
}

// Function to calculate fluid x-velocity in a single point in the pipe
double vx_pipe(point single_point_, input params_, vector<double> bessel_zeros_, vector<double> J_1_,
               vector<vector<double>> J_0_, int r_i_, int t_step_) {

    double _fb_sum = fourier_bessel(move(bessel_zeros_), move(J_1_), move(J_0_), r_i_, t_step_,
                                    params_, single_point_.getvisc());
    return ((params_.getdp() / (4 * single_point_.getvisc() * params_.getl())) * (pow(params_.getR(), 2) -
            pow(single_point_.getr(), 2))) - _fb_sum;
}
