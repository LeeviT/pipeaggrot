#ifndef ORIENTATION_AO_2012
#define ORIENTATION_AO_2012
///The orientation class -- 3d orientation model (fokker-planck)
///Direct approach:: disrectization in spherical coordinates: control volume scheme
///Many functions defined in header since template functions ...

#include <vector>
#include "boost/multi_array.hpp"
#include "triple.hpp"
#include <fstream>
#include <iostream>
#ifndef PI
#define PI  3.14159265358979323846
#endif

typedef boost::multi_array<double, 2>  d2vec;
typedef boost::detail::multi_array::multi_array_view<double, 2> d2vec_view; 

class Orientation_ao {
private:
    int s1, s2, runalready, N, M;
    ///Paramters of the system
    double k_aspect, diffusion_coeff, shear_rate;
    ///The matrices 
    d2vec W_ij, W_ipj, W_imj, W_ijp, W_ijm, area;
    std::vector<double> phi, theta;
    double h_phi, h_theta, tot_area;
    ///Set k_aspect from aspect ratio
    void set_k_aspect(const double aspect_ratio_);
    ///Get rotation velocities according to shear rate (Jeffery equation)
    ///Currently this is only used in initialize_abcde, so speed is not relevant
    ///dot_p in spherical coordinates
    std::pair<double, double> dot_p_angles(double phi_, double theta_, double shear_rate_);
    std::pair<double, double> dot_p_angles2(double phi_, double theta_);
    std::pair<double, double> dot_p_angles3(double phi_, double theta_);
public:
    // Initializer function
    void initialize(double aspect_ratio_, double shear_rate_, double diffusion_coeff_, int M_, int N_,
                    int s1_, int s2_);

    /// Initialize the matrices
    /// has to be called each time the shear rate changes (currently constant)
    /// PARAMETRIZE IF THIS HAS TO BE CALLED MORE OFTEN !!
    void initialize_abcde(double shear_rate_);

    // Update variable shear_rate
    void update_shear_rate(double shear_rate_) {
        shear_rate = shear_rate_;
    }

    ///Return an approximately uniform distribution
    d2vec get_uniform_dist();

    ///Template function here because of "array_view reference" != "array reference"
    ///Should now also eat vectors and arrays, because size information is not asked for
    template <class d2vec_view1> void get_dw_dt(d2vec_view1 &dw_, const d2vec_view1 &w_old_);

    ///since we need w_dist, this is currently a template function ..
    template <class d2vec_view1> double get_a2(const d2vec_view1 &w_dist_, int index_1_, int index_2_);
    
    ///A4 comp
    template <class d2vec_view1> double get_a4_comp(const d2vec_view1 &w_dist_);

    ///Return a4:D
    ///TODO:: Optimize this heavily
    ///since we need w_dist, this is currently a template function ..
    template <class d2vec_view1> d2vec get_tau_raw(const d2vec_view1 &w_dist_, double shear_rate_);

    template <class d2vec_view1> d2vec get_stress1(const d2vec_view1 &w_dist_);

    template <class d2vec_view1> double get_visc_part(d2vec_view1 &w_dist_);

    template <class d2vec_view1> void normalize(d2vec_view1 &w_dist_);
    
    template <class d2vec_view1> double get_tot_prob(d2vec_view1 &w_dist_);

    // Get the "amount of particles" in this distribution
    // TODO: Implement this function
    template <class d2vec_view1> double get_C(d2vec_view1 &w_dist_);

    void area_test();
    
    ///Print the current distribution 
    ///TODO take the output stream as input parameter
    template <class d2vec_view1> void shoutbox(const d2vec_view1 &w_dist_, int file_index_, int timestep);
    template <class d2vec_view1> void dist_avgs(const d2vec_view1 &w_dist_, int file_index_, int timestep);

    template <class d2vec_view1> double get_Np(const d2vec_view1 &dist_now_);

    ///Selfcheck -- mainly for Debugging
    void selfcheck();    
};
#endif
