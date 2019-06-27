#include <iostream>
#include <cstdio>
#include <cmath>
#include <fstream>
#include <array>
#include <ctime>
#include <cstdlib>
#include <vector>

// #include "main.cpp"

#include <boost/math/special_functions/bessel.hpp>
using namespace std;


// Calcute m zeros of the first kind Bessel functions of order n=1
vector<double> calc_bessel_zeros(int n) {
    vector<double> zeros(n);

    boost::math::cyl_bessel_j_zero<double>(0.0, 1, n, zeros.begin());
    return zeros;
}

// Function to calculate fluid x-velocity in a single point in the pipe
/*double vx_pipe(point single_point, double dp, double l, double radius) {
    return (dp/(4* single_point.getvisc()*l))*(pow(radius, 2) - pow(single_point.getr(), 2));
}*/

int main() {
    int n_bessel_zeros = 20;
    vector<double> bessel_zeros(n_bessel_zeros);

    bessel_zeros = calc_bessel_zeros(n_bessel_zeros);

    return 0;
}
