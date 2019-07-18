//The orientation class -- 3d orientation model (fokker-planck)
//Direct approach:: disrectization in spherical coordinates: control volume scheme
//-- working version

#include "orientation_control_vol.hpp"
#include "boost/multi_array.hpp"
#include <vector>
#include <math.h>
#include "triple.hpp"
#define PI  3.14159265358979323846

using namespace std;
using namespace boost;
using namespace boost::detail::multi_array;

namespace {
    double doublepoint_operator(const double in[][3], const triple &t1, const triple &t2) {
        double out = 0;
        out += in[0][0]*t1.x*t2.x;
        out += in[0][1]*t1.x*t2.y;
        out += in[0][2]*t1.x*t2.z;
        out += in[1][0]*t1.y*t2.x;
        out += in[1][1]*t1.y*t2.y;
        out += in[1][2]*t1.y*t2.z;
        out += in[2][0]*t1.z*t2.x;
        out += in[2][1]*t1.z*t2.y;
        out += in[2][2]*t1.z*t2.z;
        return out;
    }
}

//----------------------------------------------------------------------
/// class Orientation_ao
///PARAMETERS::
//int s1, s2;
//double  k_aspect, diffusion_coeff, shear_rate;
//Parameters related to discretization and finite difference method
//d2vec W_ij, W_ipj, W_imj, W_ijp, W_ijm, area;
//int N,M; N=phi grid size, M==theta grid size 
//std::vector<double> phi, theta;
//double h_phi, h_theta, tot_area;
//----------------------------------------------------------------------
///PRIVATE FUNCTIONS::
//Get k_aspect from aspect ratio
void Orientation_ao::set_k_aspect(const double aspect_ratio_) {
    double _aspect_two = pow(aspect_ratio_, 2);
    k_aspect = (_aspect_two - 1.0) / (_aspect_two + 1.0);
}

//dot_p in spherical coordinates 
//--do not return nans put -1e-5 instead -- careful here...
//first=dot_p_theta -- second=dot__phi
pair<double, double> Orientation_ao::dot_p_angles(double phi_, double theta_, double shear_rate_) {
    // x y z
    triple _r;
    _r.x = sin(theta_) * cos(phi_);
    _r.y = sin(theta_) * sin(phi_);
    _r.z = cos(theta_);
    // phi unit vector
    triple _p;
    _p.x = -sin(phi_);
    _p.y = cos(phi_);
    _p.z = 0;
    // theta unit vector
    triple _t;
    _t.x = cos(theta_) * cos(phi_);
    _t.y = cos(theta_) * sin(phi_);
    _t.z = -sin(theta_);
    pair<double, double> _dot_p_s;
    double _kappa[3][3] = {{0.0, 0.0, 0.0,}, {0.0,0.0,0.0}, {0.0,0.0,0.0}};
    double _kappa_t[3][3] = {{0.0, 0.0, 0.0,}, {0.0,0.0,0.0}, {0.0,0.0,0.0}};
    _kappa[s2][s1] = _kappa_t[s1][s2] = shear_rate_;
    double _lambda = k_aspect;
    double _lambda_p = (_lambda + 1) / 2.0;
    // now lambda m == 0
    double _lambda_m = (_lambda - 1) / 2.0;
    double _k_rp = doublepoint_operator(_kappa, _r, _p);
    double _kt_rp = doublepoint_operator(_kappa_t, _r, _p);
    double _k_rt = doublepoint_operator(_kappa, _r, _t);
    double _kt_rt = doublepoint_operator(_kappa_t, _r, _t);
    _dot_p_s.second  = _lambda_m * _kt_rp + _lambda_p * _k_rp;
    _dot_p_s.first  = _lambda_m * _kt_rt + _lambda_p * _k_rt;
    return _dot_p_s;
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------
//----------------------------------------------------------------------
// PUBLIC FUNCTIONS::
// initialize these --
// has to be called each time the shear rate changes (currently constant)
// PARAMETRIZE IF THIS HAS TO BE CALLED MORE OFTEN !!
// Since Only thing that actually changes 
// upon shear rate modification is  dot_p_angles
void Orientation_ao::initialize_abcde(double shear_rate_) {
    // in W_???: "p" refers to +1 and "m" to i-1  (imj->i-1,j) etc.
    W_ij.resize(extents[N][M]);
    W_ijp.resize(extents[N][M]);
    W_ijm.resize(extents[N][M]);
    W_ipj.resize(extents[N][M]);
    W_imj.resize(extents[N][M]);
    area.resize(extents[N][M]);
    tot_area = 0.0;
    double _phi, _phi_m, _phi_p, _theta, _theta_m, _theta_p;
    for (int i = 0; i < N; i++) {
        _phi = phi[i];
        _phi_p = phi[i] + 0.5*h_phi;
        _phi_m = phi[i] - 0.5*h_phi;
        for (int j = 0; j < M; j++) {
            _theta = theta[j];
            _theta_p = theta[j] + 0.5*h_theta;
            _theta_m = theta[j] - 0.5*h_theta;
            // Area is exactly given by this formulation:
            area[i][j] = h_phi * fabs(cos(_theta_m) - cos(_theta_p));
            if (_theta_m < PI/2 && _theta_p > PI/2)
                area[i][j] = h_phi * (2 - (1 - fabs(cos(_theta_m))) - (1 - fabs(cos(_theta_p))));
                // Orientation part (same below)
                W_ipj[i][j] = (-sin(_theta) / 2 * dot_p_angles(_phi_p, _theta, shear_rate_).second / sin(_theta)
                // Diffusion/aggregation part (same below)
                + diffusion_coeff / (sin(_theta) * h_phi))
                * (h_theta / area[i][j]);
                W_imj[i][j] = (sin(_theta) / 2 * dot_p_angles(_phi_m, _theta, shear_rate_).second / sin(_theta)
                + diffusion_coeff / (sin(_theta) * h_phi))
                * (h_theta / area[i][j]);
            W_ijp[i][j] = (-sin(_theta_p) / 2 * dot_p_angles(_phi, _theta, shear_rate_).first
                + diffusion_coeff * sin(_theta_p) / h_theta )
                * (h_phi / area[i][j]);
            W_ijm[i][j] = (sin(_theta_m) / 2 * dot_p_angles(_phi, _theta_m, shear_rate).first
                + diffusion_coeff * sin(_theta_m) / h_theta)
                * (h_phi / area[i][j]);
            W_ij[i][j] = W_ijp[i][j] +  W_ijm[i][j] + W_ipj[i][j] + W_imj[i][j]
                - 2 * diffusion_coeff * sin(_theta_m) * h_phi / (h_theta * area[i][j])
                - 2 * diffusion_coeff * sin(_theta_p) * h_phi / (h_theta * area[i][j])
                - 4 * diffusion_coeff * h_theta / (sin(_theta) * h_phi * area[i][j]);
            tot_area += area[i][j];
        }
    }
}

void Orientation_ao::initialize(const double aspect_ratio_, double shear_rate_, double diffusion_coeff_,
                                int N_, int M_, int s1_, int s2_) {
    runalready = 0; // this is for deciding whether to create new files or append
    s1 = s1_;
    s2 = s2_;
    M = M_;
    N = N_; // system size
    h_phi = 2 * PI / N;
    phi.resize(N);
    h_theta = PI/(M);
    theta.resize(M);
    for (int i = 0; i < N; i++) {
        phi[i] = i * h_phi + 0.5 * h_phi;
    }
    for (int i = 0; i < M; i++) {
        theta[i] = h_theta * i + 0.5 * h_theta;
    }
    shear_rate = shear_rate_; // shear rate
    set_k_aspect(aspect_ratio_); // aspect ratio initialization
    diffusion_coeff = diffusion_coeff_; // diffusion coefficient
    initialize_abcde(shear_rate_);
}
// Return an approximately uniform distribution
d2vec Orientation_ao::get_uniform_dist() {
    d2vec _w_uniform(extents[N][M]);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            _w_uniform[i][j] = 1.0 / (tot_area);
        }
    }
    return _w_uniform;
}
// Get dw_dt --- in header...
// Template function here because of "array_view reference" != "array reference"
// Functions defined in header ...
// Get dw_dt
// Template function here because of "array_view reference" != "array reference"
// Should now also eat vectors and arrays, because size information is not asked for
template <class d2vec_view1> void Orientation_ao::get_dw_dt(d2vec_view1 &dw_, const d2vec_view1 &w_old_) {
    int _im1, _ip1, _jm1, _jp1;
    for (int i = 0; i < N; i++) {
        _im1 = i - 1; _ip1 = i + 1; if(!i) _im1 = N - 1; else if(i == N - 1) _ip1 = 0;
        for (int j = 1; j < M - 1; j++) {
            _jm1 = j - 1; _jp1 = j + 1;
            dw_[i][j] += w_old_[i][j] * W_ij[i][j] + w_old_[i][_jp1] * W_ijp[i][j]
                       + w_old_[i][_jm1] * W_ijm[i][j] + w_old_[_ip1][j] * W_ipj[i][j] + w_old_[_im1][j] * W_imj[i][j];
        }
        // FOR control volume scheme:: flux over the poles is zero...
        // SPECIAL CASE:: NORTH POLE
        dw_[i][0] += w_old_[i][0] * W_ij[i][0] + w_old_[i][1] * W_ijp[i][0]
                   + w_old_[_ip1][0] * W_ipj[i][0] + w_old_[_im1][0] * W_imj[i][0];
        // SPECIAL CASE:: SOUTH POLE
        dw_[i][M - 1] += w_old_[i][M - 1] * W_ij[i][M - 1] + w_old_[i][M - 2] * W_ijm[i][M - 1]
                       + w_old_[_ip1][M - 1] * W_ipj[i][M - 1] + w_old_[_im1][M - 1] * W_imj[i][M - 1];
    }
}

// GET a2
// since we need w_dist, this is currently a template function...
template <class d2vec_view1> double Orientation_ao::get_a2(const d2vec_view1 &w_dist_, int index_1_, int index_2_) {
    // Get "a2" second order orientation tensor (function of w_dist)
    double _a2[3][3] = {0.0};
    for (int phi_now = 0; phi_now < N; phi_now++) {
        for (int theta_now = 0; theta_now < M; theta_now++) {
            triple _p_triple(phi[phi_now], theta[theta_now]);  // These should be stored somewhere
            double _p[3] = {_p_triple.x, _p_triple.y, _p_triple.z}; // or even these
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    // basically one could even store something like a4_raw[phi_now][theta_now][i][j][k][l]
                    _a2[i][j] += _p[i] * _p[j] * area[phi_now][theta_now] * w_dist_[phi_now][theta_now];
                }
            }
        }
    }
    return _a2[index_1_][index_2_];
}

// A4 comp
template <class d2vec_view1> double Orientation_ao::get_a4_comp(const d2vec_view1 &w_dist_) {
    // Get "a2" second order orientation tensor (function of w_dist)
    double _a2[3][3] = {0.0};
    for (int phi_now = 0; phi_now < N; phi_now++) {
        for (int theta_now = 0; theta_now < M; theta_now++) {
            triple _p_triple(phi[phi_now], theta[theta_now]);  // These should be stored somewhere
            double _p[3] = {_p_triple.x, _p_triple.y, _p_triple.z}; // or even these
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    // basically one could even store something like a4_raw[phi_now][theta_now][i][j][k][l]
                    _a2[i][j] += _p[i] * _p[j] * area[phi_now][theta_now] * w_dist_[phi_now][theta_now];
                }
            }
        }
    }
    double _a4_comp[3][3][3][3] = {0.0};
    double _a4[3][3][3][3] = {0.0};
    for (int phi_now = 0; phi_now < N; phi_now++) {
        for (int theta_now = 0; theta_now < M; theta_now++) {
            triple _p_triple(phi[phi_now], theta[theta_now]);  // These should be stored somewhere
            double _p[3] = {_p_triple.x, _p_triple.y, _p_triple.z}; // or even these
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    for (int k = 0; k < 3; k++) {
                        for (int l = 0; l < 3; l++) {
                            // basically one could even store something like a4_raw[phi_now][theta_now][i][j][k][l]
                            _a4[i][j][k][l] += _p[i] * _p[j] * _p[k] * _p[l] * area[phi_now][theta_now]
                                             * w_dist_[phi_now][theta_now];
                        }
                    }
                }
            }
        }
    }
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                for (int l = 0; l < 3; l++) {
                    _a4_comp[i][j][k][l] = _a2[i][j] * _a2[k][l];
                }
            }
        }
    }
    std::cout << "a4 \t a4_comp" << "\n";
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                for (int l = 0; l < 3; l++) {
                    std::cout << _a4[i][j][k][l] << "\t" << _a4_comp[i][j][k][l] <<
                              "\t" << _a4[i][j][k][l] - _a4_comp[i][j][k][l] << "\n";
                }
            }
        }
    }
    std::cout << std::endl;
    exit(1);
}

// Return a4:D
// TODO: Optimize this heavily
// since we need w_dist, this is currently a template function ..
template <class d2vec_view1> d2vec Orientation_ao::get_tau_raw(const d2vec_view1 &w_dist_, double shear_rate_) {
    d2vec _a4_D(boost::extents[3][3]);
    std::fill(_a4_D.data(), _a4_D.data() + _a4_D.num_elements(), 0.0);
    // D should be defined somewhere else
    double _D[3][3] = {0.0};
    _D[s1][s2] = _D[s2][s1] = shear_rate_ / 2.0;
    // Get "a4" fourth order orientation tensor (function of w_dist)
    double _a4[3][3][3][3] = {0.0};
    for (int phi_now = 0; phi_now < N; phi_now++) {
        for (int theta_now = 0; theta_now < M; theta_now++) {
            triple _p_triple(phi[phi_now], theta[theta_now]); // These should be stored somewhere
            double _p[3] = {_p_triple.x, _p_triple.y, _p_triple.z}; // or even these
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    for (int k = 0; k < 3; ++k) {
                        for (int l = 0; l < 3; ++l) {
                            // basically one could even store something like a4_raw[phi_now][theta_now][i][j][k][l]
                            _a4[i][j][k][l] += _p[i] * _p[j] * _p[k] * _p[l] * area[phi_now][theta_now]
                                             * w_dist_[phi_now][theta_now];
                        }
                    }
                }
            }
        }
    }
    // Get a4:D
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                for (int l = 0; l < 3; l++) {
                    _a4_D[i][j] += _a4[i][j][k][l] * _D[k][l];
                }
            }
        }
    }
    return _a4_D;
}

template <class d2vec_view1> d2vec Orientation_ao::get_stress1(const d2vec_view1 &w_dist_) {
    std::cout << "GET STRESS " << std::endl;
    double _Np = 300; // TODO: put this somewhere else (idk??)
    // D should be defined somewhere else
    double _D[3][3] = {0.0};
    _D[s1][s2] = _D[s2][s1] = shear_rate / 2.0;
    d2vec _a4_D = get_tau_raw(w_dist_, shear_rate);
    d2vec _stress(boost::extents[3][3]);

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            _stress[i][j] = 2 * (_D[i][j] + _Np * _a4_D[i][j]);
            std::cout << "#" << i << j << " " << _stress[i][j] << std::endl;
        }
    }
    return _stress;
}

template <class d2vec_view1> double Orientation_ao::get_visc_part(d2vec_view1 &w_dist_) {
    normalize(w_dist_); // w_dist_ is not const because of this line...
    double _viscosity_part, _a4[3][3][3][3] = {0.0};
    for (int phi_now = 0; phi_now < N; phi_now++) {
        for (int theta_now = 0; theta_now < M; theta_now++) {
            triple _p_triple(phi[phi_now], theta[theta_now]);  // These should be stored somewhere
            double _p[3] = {_p_triple.x, _p_triple.y, _p_triple.z}; // or even these
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    for (int k = 0; k < 3; k++) {
                        for (int l = 0; l < 3; l++) {
                            _a4[i][j][k][l] += _p[i] * _p[j] * _p[k] * _p[l] * area[phi_now][theta_now]
                                              * w_dist_[phi_now][theta_now];
                        }
                    }
                }
            }
        }
    }
    _viscosity_part = 2 * _a4[s1][s2][s1][s2];
    return _viscosity_part;
}

template <class d2vec_view1> void Orientation_ao::normalize(d2vec_view1 &w_dist_) {
    double _tot_prob = 0.0;
    for (int phi_now = 0; phi_now < N; phi_now++) {
        for (int theta_now = 0; theta_now < M; theta_now++) {
            _tot_prob += w_dist_[phi_now][theta_now] * area[phi_now][theta_now];
        }
    }
    std::cout << "#tot_prob1= " << _tot_prob << std::endl;
}

template <class d2vec_view1> double Orientation_ao::get_tot_prob(d2vec_view1 &w_dist_) {
    double _tot_prob = 0.0;
    for (int phi_now = 0; phi_now < N; phi_now++) {
        for (int theta_now = 0; theta_now < M; theta_now++) {
            _tot_prob += w_dist_[phi_now][theta_now] * area[phi_now][theta_now];
        }
    }
    return _tot_prob;
}

// Get the "amount of particles" in this distribution
template <class d2vec_view1> double Orientation_ao::get_C(d2vec_view1 &w_dist_) {
    // TODO: Implement this function
}

void Orientation_ao::area_test() {
    double _tot_area = 0.0;
    for (int phi_now = 0; phi_now < N; phi_now++) {
        for (int theta_now = 0; theta_now < M; theta_now++) {
            _tot_area += area[phi_now][theta_now];
        }
    }
    //~ std::cerr << tot_area/(4*PI) << std::endl;
}

// Print the current distribution
// TODO: take the output stream as input parameter
template <class d2vec_view1> void Orientation_ao::shoutbox(const d2vec_view1 &w_dist_, int file_index_, int timestep_) {
    std::stringstream _filename;
    // _filename << "class_" << file_index_ << "_t_" << timestep_ << ".dat";
    _filename << "class_" << file_index_ << ".dat";
    std::ofstream _outfile;
    if (!runalready) _outfile.open(_filename.str().c_str());
    else _outfile.open(_filename.str().c_str(), std::ios::out | std::ios::app);
    runalready++;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            triple _tmp(phi[i], theta[j]);
            _tmp.add_r(w_dist_[i][j] / (4 * PI));
            // _tmp.shoutbox_file(_outfile);
        }
        _outfile << "\n";
    }
    // Print first again in order to get grid in gnuplot
    for (int j = 0; j < M; j++) {
        triple _tmp(phi[0], theta[j]);
        _tmp.add_r(w_dist_[0][j] / (4 * PI));
        _tmp.shoutbox_file(_outfile);
    }
    _outfile << "\n";
    _outfile << std::endl;
}

// Print average theta and phi values of the distribution
template <class d2vec_view1> void Orientation_ao::dist_avgs(const d2vec_view1 &w_dist_, int file_index_,
                                                            int timestep_) {
    std::stringstream _filename;
    // _filename << "dist/class_" << file_index_ << "_" << timestep_ << ".dat";
    _filename << "dist/class_" << file_index_ << ".dat";
    std::ofstream _outfile;
    // _outfile.open(_filename.str().c_str());
    _outfile.open(_filename.str().c_str(), std::ios::out | std::ios::app);
    runalready++;
    double phi_sum = 0.0, theta_sum = 0.0;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            // _outfile << w_dist_[i][j] / (4 * PI) << "\t" << phi[i] << "\t" << theta[j] << std::endl;
            // _outfile << w_dist_[i][j] * phi[i] * (4 * PI) << "\t" << w_dist_[i][j] * theta[j] * (4 * PI) << endl;
            phi_sum += w_dist_[i][j] * phi[i] * (4 * PI);
            theta_sum += w_dist_[i][j] * theta[i] * (4 * PI);
        }
    }
    _outfile << timestep_ << "\t" << phi_sum / (20.0 * 20.0) << "\t" << theta_sum / (20.0 * 20.0) << endl;
}

template <class d2vec_view1> double Orientation_ao::get_Np(const d2vec_view1 &dist_now_) {
    double _Np_now = 0.0;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            _Np_now += dist_now_[i][j] * area[i][j];
        }
    }
    return _Np_now;
}

// Selfcheck -- mainly for Debugging
void Orientation_ao::selfcheck() {
    double _dt_unity;

    cout << "\n\n";
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            cout << (phi[i] / (2 * PI)) * 360 << " " << (theta[j] / (2 * PI)) * 360 << " "
                 << dot_p_angles(phi[i], theta[j], shear_rate).second << endl;
        }
    }
    cout << "\n\n";
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            cout << (phi[i] / (2 * PI)) * 360 << " " << (theta[j] / (2 * PI)) * 360 << " "
                 << dot_p_angles(phi[i], theta[j], shear_rate).first << endl;
        }
    }
    cout << "\n\n";
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            cout << (phi[i] / (2 * PI)) * 360 << " " << (theta[j] / (2 * PI)) * 360 << " " << W_ij[i][j] << endl;
        }
    }
    cout << "\n\n";
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            cout << (phi[i] / (2 * PI)) * 360 << " " << (theta[j] / (2 * PI)) * 360 << " " << W_ijp[i][j] << endl;
        }
    }
    cout << "\n\n";
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            cout << (phi[i] / (2 * PI)) * 360 << " " << (theta[j] / (2 * PI)) * 360 << " " << W_ijm[i][j] << endl;
        }
    }
    cout << "\n\n";
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            cout << (phi[i] / (2 * PI)) * 360 << " " << (theta[j] / (2 * PI)) * 360 << " " << W_ipj[i][j] << endl;
        }
    }
    cout << "\n\n";
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            cout << (phi[i] / (2 * PI)) * 360 << " " << (theta[j] / (2 * PI)) * 360 << " " << W_imj[i][j] << endl;
        }
    }
    cout << "\n\n";
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            _dt_unity = W_ij[i][j] + W_ijp[i][j] + W_ijm[i][j] + W_ipj[i][j] + W_imj[i][j];
            cout << (phi[i] / (2 * PI)) * 360 << " " << (theta[j] / (2 * PI)) * 360 << " " << _dt_unity << endl;
        }
    }
    cout << "\n\n";
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            triple _temp(phi[i], theta[j]);
            _dt_unity = W_ij[i][j] + W_ijp[i][j] + W_ijm[i][j]; //+W_ipj[i][j]+ W_imj[i][j];
            cout << _temp.x << " " << _temp.y << " " << _temp.z << " " << _dt_unity << endl;
        }
    }
    cout << "\n\n";

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            _dt_unity = W_ij[i][j] + W_ijp[i][j] + W_ijm[i][j] + W_ipj[i][j] + W_imj[i][j];
            if (j == 0) {
                triple _temp1(phi[i], theta[j]);
                _temp1.add_r(_dt_unity);
                cout << _temp1.x << "\t" << _temp1.y << "\t" << _temp1.z << endl;
            } else if (j == M) {
                triple _temp2(phi[i], theta[j]);
                _temp2.add_r(_dt_unity);
                cout << _temp2.x << "\t" << _temp2.y << "\t" << _temp2.z << endl;
            } else {
                triple _temp(phi[i], theta[j]);
                _temp.add_r(_dt_unity);
                _temp.transform();
                cout << _temp.x << "\t" << _temp.y << "\t" << _temp.z << endl;
            }
        }
    }
    cout << "\n\n";
    
    // Check dot_p_angles -- 2 cases
    cout << "\n\n";
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            triple _temp(phi[i], theta[j]);
            pair<double, double> _temp_dt_angles = dot_p_angles(phi[i], theta[j], shear_rate);
            _temp.add_r(_temp_dt_angles.first);
            cout << _temp.x << "\t" << _temp.y << "\t"<< _temp.z << endl;
        }             
    }
    cout << "\n\n";
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            triple _temp(phi[i], theta[j]);
            pair<double, double> _temp_dt_angles = dot_p_angles(phi[i], theta[j], shear_rate);
            _temp.add_r(_temp_dt_angles.second);
            cout << _temp.x << "\t" << _temp.y << "\t"<< _temp.z << endl;
        }             
    }
    cout << "\n\n";
    // Check dot_p recreate -- 3 cases
    cout << "\n\n";
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            triple _temp(phi[i], theta[j]);
            pair<double, double> _temp_dt_angles = dot_p_angles(phi[i], theta[j], shear_rate);
            double _px = cos(theta[j]) * cos(phi[i]) * _temp_dt_angles.first
            - sin(theta[j]) * sin(phi[i]) * _temp_dt_angles.second;
            _temp.add_r(_px);
            cout << _temp.x << "\t" << _temp.y << "\t"<< _temp.z << endl;
        }             
    }
    cout << "\n\n";
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            triple _temp(phi[i], theta[j]);
            pair<double, double> _temp_dt_angles = dot_p_angles(phi[i], theta[j], shear_rate);
            double _py = cos(theta[j]) * sin(phi[i]) * _temp_dt_angles.first
            + sin(theta[j]) * cos(phi[i]) * _temp_dt_angles.second;
            _temp.add_r(_py);
            cout << _temp.x << "\t" << _temp.y << "\t"<< _temp.z << endl;
        }             
    }
    cout << "\n\n";
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            triple _temp(phi[i], theta[j]);
            pair<double, double> _temp_dt_angles = dot_p_angles(phi[i], theta[j], shear_rate);
            double _pz = -sin(theta[j]) * _temp_dt_angles.first;
            _temp.add_r(_pz);
            cout << _temp.x << "\t" << _temp.y << "\t"<< _temp.z << endl;
        }             
    }
    cout << "\n\n";
}

template void Orientation_ao::get_dw_dt<boost::detail::multi_array::multi_array_view<double, 2>>
              (d2vec_view&_, const d2vec_view&);
template void Orientation_ao::normalize<boost::detail::multi_array::multi_array_view<double, 2>>
              (d2vec_view&_);
template void Orientation_ao::shoutbox<boost::detail::multi_array::multi_array_view<double, 2>>
              (const d2vec_view&_, int file_index_, int timestep_);
template void Orientation_ao::dist_avgs<boost::detail::multi_array::multi_array_view<double, 2>>
        (const d2vec_view&_, int file_index_, int timestep_);
template double Orientation_ao::get_a2<boost::detail::multi_array::multi_array_view<double, 2>>
                (const d2vec_view&_, int index_1_, int index_2_);
template double Orientation_ao::get_a4_comp<boost::detail::multi_array::multi_array_view<double, 2>>
                (const d2vec_view&_);
template double Orientation_ao::get_visc_part<boost::detail::multi_array::multi_array_view<double, 2>>
                (d2vec_view&_);
template double Orientation_ao::get_tot_prob<boost::detail::multi_array::multi_array_view<double, 2>>
                (d2vec_view&_);
template double Orientation_ao::get_C<boost::detail::multi_array::multi_array_view<double, 2>>
                (d2vec_view&_);
template double Orientation_ao::get_Np<boost::detail::multi_array::multi_array_view<double, 2>>
                (const d2vec_view&_);
template d2vec Orientation_ao::get_tau_raw<boost::detail::multi_array::multi_array_view<double, 2>>
               (const d2vec_view&_, double shear_rate_);
template d2vec Orientation_ao::get_stress1<boost::detail::multi_array::multi_array_view<double, 2>>
               (const d2vec_view&_);