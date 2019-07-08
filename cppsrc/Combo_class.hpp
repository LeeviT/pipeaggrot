#ifndef COMBOCLASS_AO_2012
#define COMBOCLASS_AO_2012

#include "orientation_control_vol.hpp"
// #include "orientation_combo.hpp"
#include "dist_combo.hpp"
#include "point.hpp"
#include "io.hpp"
#include <vector>

typedef boost::multi_array<double, 3>  d3vec;

class Combo_class {
    double shear_rate, diffusion_coeff_max, visco0, Np_max, beta;
    int classes, M, N;
    double d_classes; 
    d3vec distribution;
    d3vec dist_dt;
    Dist_ao aggregation_class;
    std::vector<Orientation_ao> orientation_class;
    public:
    ///INITIALIZE
    void initialize(double shear_rate_, input params);
    ///GET DT
    void get_dt(double* y_dot, double* y_now);
    ///SET INITIAL UNIFORM STATE
    void set_uniform_state(double* y_now);
    ///LOAD STATE
    void load_state(double* y_now);
    ///SAVE THE CURRENT STATE
    void save_state();
    ///Get visc raw
    double get_visc_raw(int i1_, int i2_);
    ///Get tau raw
    double get_tau_raw(int i1_, int i2_);
    ///Get tot prob
    void get_tot_prob(double time_now = -1.0);
	void save_aggr_distribution(double time_now = -1.0);
    void update_shear_rate(double shear_rate);

    //EVERYTHING BELOW NEEDS CHECKING
/*
    ///NOT NEEDED YET
    ///PRINT CURRENT STATE --- should be improved
    void print_out();
    ///SHOUTBOX PRINTOUT MODE
    void shoutbox();
    ///GET visc part
    void get_visc_part(double time);
    ///GET visc full
    double get_visc_tot(int s1_, int s2_);
    ///GET A2 part
    void get_a2_(int index1, int index2, double time_);
    ///Get the size distribution (piecewise for each class)
    void get_C();
    //DEBUGING...
    void get_a4_comp();
*/
};

#endif
