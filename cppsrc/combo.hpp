#ifndef COMBOCLASS_AO_2012
#define COMBOCLASS_AO_2012

#include "orientation.hpp"
#include "distribution.hpp"
#include "point.hpp"
#include "io.hpp"

#include <vector>

typedef boost::multi_array<double, 3>  d3vec;

class combo {
    double shear_rate, diffusion_coeff_max, Np_max, visco0, beta;
    int classes, M, N;
    double d_classes; 
    d3vec distribution;
    d3vec dist_dt;
    Dist_ao aggregation_class;
    std::vector<Orientation_ao> orientation_class;
    public:
    ///INITIALIZE
    void initialize(double shear_rate_, input params_);
    ///GET DT
    void get_dt(double * y_dot_, double * y_now_);
    ///SET INITIAL UNIFORM STATE
    void set_uniform_state(double * y_now_);
    ///LOAD STATE
    void load_state(double * y_now_);
    ///SAVE THE CURRENT STATE
    void save_state();
    ///Get visc raw
    double get_visc_raw(int i1_, int i2_, int timestep_);
    ///Get tau raw
    double get_tau_raw(int i1_, int i2_);
    ///Get tot prob
    void get_tot_prob(double time_now_ = -1.0);
	void save_aggr_distribution(double time_now_ = -1.0);
    void update_shear_rate(double shear_rate_);
};

#endif
