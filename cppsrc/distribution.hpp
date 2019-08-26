#ifndef DIST_AO_2012
#define DIST_AO_2012

#include "boost/multi_array.hpp"

typedef boost::multi_array<double, 1> d1vec;

// Class to handle distribution stuff
class Dist_ao {
private:
    int system_size;
    double dr, tau, visc0, v_d_max, v_c_max, shear_rate_max, v_d, v_c, beta_now;

    template <class d1vec_view1> double w_dist_dt_aggregation(int n_, const d1vec_view1 &w_dist_) {
        //Assuming constant spacing for dr!!
        double _w_dist_dt = 0;
        for (int i = 0; i < n_; i++) {
            _w_dist_dt += v_c * w_dist_[i] * dr;
        }
        for (int i = n_ + 1; i < system_size; i++) {
            _w_dist_dt += v_d * w_dist_[i] * dr;
        }
        for (int i = 0; i < n_; i++) {
            _w_dist_dt -= v_d * w_dist_[n_] * dr;
        }
        for (int i = n_ + 1; i < system_size; i++) {
            _w_dist_dt -= v_c * w_dist_[n_] * dr;
        }
        return _w_dist_dt;
    }

public:
    // Initialize class instance
    void initialize(int wdist_size_, double shear_rate_, double beta_, double shear_rate_max_) {
        beta_now = beta_;
        system_size = wdist_size_;
        shear_rate_max = shear_rate_max_;
        v_c_max = 1.0;
        v_d_max = v_c_max * beta_now;
        system_size = wdist_size_;
        dr = 1.0 / wdist_size_;
        update_v_aggregation(shear_rate_);
    }
    // Get w_dist_dt (only aggregation and fragmentation)
    template <class d1vec_view1> void get_wdist_dt(d1vec_view1 &w_dist_dt_, const d1vec_view1 &w_dist_) {
        for (int i = 0; i < system_size; i++) {
            w_dist_dt_[i] += w_dist_dt_aggregation(i, w_dist_);
        }
    }
    // Update aggregation and breakage velocities according to shear rate
    void update_v_aggregation(const double shear_rate_) {
        if (shear_rate_ > shear_rate_max) {
            v_d = v_d_max;
            v_c = 0;
        } else {
            v_d = v_d_max * (shear_rate_ / shear_rate_max);
            v_c = v_c_max * (shear_rate_max - shear_rate_) / shear_rate_max;
        }
    }
    // Getter for beta
    double get_beta() {
        return beta_now;
    }
};
#endif
