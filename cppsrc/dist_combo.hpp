#ifndef DIST_AO_2012
#define DIST_AO_2012
///The distribution class -- aggregation orientation model
#include "boost/multi_array.hpp"
typedef boost::multi_array<double, 1> d1vec;

class Dist_ao {
private:
    int system_size;
    double dr, tau, visc0;
    double v_d_max, v_c_max, shear_rate_max;
    double v_d, v_c;
    double beta_now;
    ///update aggregation and breakage velocities according to shear rate
    void update_v_aggregation(const double shear_rate) {
        if (shear_rate > shear_rate_max) {
           v_d = v_d_max;
           v_c = 0;
        } else {
           v_d = v_d_max * (shear_rate/shear_rate_max);
           v_c = v_c_max * (shear_rate_max - shear_rate)/shear_rate_max;
        }
    }
    template <class d1vec_view1>
    double w_dist_dt_aggregation(int n, const d1vec_view1 &w_dist) {
        //Assuming constant spacing for dr!!
        double w_dist_dt = 0;
        for (int i = 0; i < n; i++) {
            w_dist_dt += v_c * w_dist[i] * dr;
        }
        for (int i = n + 1; i < system_size; i++) {
            w_dist_dt += v_d * w_dist[i] * dr;
        }
        for (int i = 0; i < n; i++) {
            w_dist_dt -= v_d * w_dist[n] * dr;
        }
        for (int i = n + 1; i < system_size; i++) {
            w_dist_dt -= v_c * w_dist[n] * dr;
        }
        return w_dist_dt;
    }
public:
    void update_shear_rate(double shear_rate) {
        update_v_aggregation(shear_rate);
    }
    ///INITIALIZE
    void initialize(int wdist_size, double shear_rate, double beta, double shear_rate_max_) {
        beta_now = beta;
        system_size = wdist_size;
        shear_rate_max = shear_rate_max_;
        v_c_max = 1.0;
        v_d_max = v_c_max * beta_now;
        system_size = wdist_size;
        dr = 1.0 / wdist_size;
        update_shear_rate(shear_rate);
    }
    ///get w_dist_dt (only aggregation and fragmentation)
    template <class d1vec_view1>
    void get_wdist_dt(d1vec_view1 &w_dist_dt, const d1vec_view1 &w_dist) {
        for (int i = 0; i < system_size; i++) {
            w_dist_dt[i] += w_dist_dt_aggregation(i, w_dist);
        }
    }
    //~ template <class d1vec_view1>
    //~ void print_dist(const d1vec_view1 &w_dist, ofstream &outfile){
        //~ for(int i=0; i!= system_size; ++i){
            //~ outfile << i << " " << w_dist[i]; 
        //~ }
    //~ }
    double get_beta() {
        return beta_now;
    }
};
#endif
