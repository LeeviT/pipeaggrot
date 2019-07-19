#include "Combo_class.hpp"
#include <algorithm>
#include <iostream>
#include <string.h> // memcpy is in string.h ...

typedef boost::multi_array_types::index_range range;
typedef d3vec::index indexb;
using namespace std;
using namespace boost;

// Now we need extra space os size "classes" for the aggregation distribution
void Combo_class::initialize(double shear_rate_, input params_) {

    double _diffusion_coeff;
    shear_rate = shear_rate_;
    diffusion_coeff_max = params_.getdiff_coeff();
    classes = params_.getclasses(); M = params_.getphi_grid(); N = params_.gettheta_grid();
    // Extra space of size classes*M for aggregation part <-- some overhead here...
    distribution.resize(extents[classes][N + 1][M]);
    dist_dt.resize(extents[classes][N + 1][M]);
    orientation_class.resize(classes); // This many orientation classes
    d_classes = 1.0 / classes;
    visco0 = params_.getvisc0();
    Np_max = params_.getNp_max();
    // By default aggregation model is not used but if condition is true, create aggregation class
    if (classes > 1)  aggregation_class.initialize(classes, shear_rate, params_.getbeta(), params_.getshear_rate_max());
    // Each size has its own orienation object
    // The diffusion coefficient is a function of the class
    if (M != 1 && N != 1) { // If this is false forget about the orientation part
        for (indexb i = 0; i < classes; i++) {
	        _diffusion_coeff = diffusion_coeff_max * (i + 1.0) / (classes * 1.0);
            orientation_class[i].initialize(params_.getaspect_ratio(), shear_rate, _diffusion_coeff, N, M,
                                            params_.gets1(), params_.gets2());
        }
	}
    if (params_.getwhat_todo() == -1) {
        for (indexb i = 0; i < classes; i++) {
            orientation_class[i].selfcheck();
        }
        throw("selfcheck");
    }
}

// After this function the distribution stored in Combo_class is the one given as input * y_now_
void Combo_class::get_dt(double * y_dot_, double * y_now_) {
    // Quite some overhead here since have to copy a large array at each step ... 
    // Working with boost multiarray -- that's why
    memcpy(distribution.data(), y_now_, sizeof(double) * distribution.num_elements());
    // Initialize dist_dt with zeroes
    fill(dist_dt.data(),dist_dt.data() + dist_dt.num_elements(), 0.0);
    // The tempView:s are created since need non-const reference to array view
    // Thus the array views cannot be created upon function call 
    // Probably there is an easier/cleaner way...
    // Orientation part 
    if (M != 1 && N != 1) { //If this is flase forget about the orientation part
        for (int i = 0; i < classes; i++) {
            multi_array_ref<double, 2>::array_view<2>::type _tempView_dt = dist_dt[indices[i][range(0, N)][range()]];
            orientation_class[i].get_dw_dt(_tempView_dt, distribution[indices[i][range(0, N)][range()]]);
        }
	}
    // Aggregation_fragmentation part
    if (classes > 1) {
	    multi_array_ref<double, 1>::array_view<1>::type _tempView_dt = dist_dt[indices[range()][N][0]];
        aggregation_class.get_wdist_dt(_tempView_dt, distribution[indices[range()][N][0]]);
    }
    // Copy the calculated dist_dt (multiarray) to y_dot_ ...
    memcpy(y_dot_, dist_dt.data(), sizeof(double) * dist_dt.num_elements());
}

// Simple function for setting uniform initial state
void Combo_class::set_uniform_state(double * y_now_) {
    // Step1:: the distribution
    fill(distribution.data(), distribution.data() + distribution.num_elements(), 0.0);
    for (int i = 0; i < classes; i++) {
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < M; k++) {
                distribution[i][j][k] = 1.0 / (4.0 * PI);
            }
        }
		distribution[i][N][0] = 1.0 / classes;
    }
    // Step2:: copy the data from distribution to "the N_Vector"
    memcpy(y_now_, distribution.data(), sizeof(double) * distribution.num_elements());
}

// Simple function for loading a state from file
void Combo_class::load_state(double * y_now_) {
    // Step1:: read the State from a file into distribution
    ifstream _read_file("saved_state.dat");
    for (int i = 0; i < classes; i++) {
        for (int j = 0; j <= N; j++) {
            for (int k = 0; k < M; k++) {
                _read_file >> distribution[i][j][k];
            }
        }
    }
    // Step2:: copy the data from distribution to "the N_Vector"
    memcpy(y_now_, distribution.data(), sizeof(double) * distribution.num_elements());
}
    
// Simple function for saving a state from a file
void Combo_class::save_state() {
    ofstream _save_file("saved_state.dat");
    for (int i = 0; i < classes; i++) {
        for (int j = 0; j <= N; j++) {
            for (int k = 0; k < M; k++) {
                _save_file << distribution[i][j][k] << " ";
            }
        }
        _save_file << "\n";
    }
}

double Combo_class::get_visc_raw(int s1_, int s2_, int timestep_) {
    for (int i = 0; i < classes; i++) {
        multi_array_ref<double, 2>::array_view<2>::type _tempView_wdist = distribution[indices[i][range(0, N)][range()]];
        orientation_class[i].dist_avgs(_tempView_wdist, i, timestep_);
    }
    return get_tau_raw(s1_, s2_) / shear_rate;
}

// This should give the stress part due to the orientation...
double Combo_class::get_tau_raw(int s1_, int s2_) {
    double _Npi, _aggr_part;

    d2vec _tau_raw(boost::extents[3][3]);

    std::fill(_tau_raw.data(), _tau_raw.data() + _tau_raw.num_elements(), 0.0); // Fill with zeroes

    for (int i = 0; i < classes; i++) {
        // Np varies linearly with the population variable n
        if (classes == 1) {
            _Npi = Np_max;
        } else {
            _Npi = Np_max * (i + 1) * 1.0 / classes;
        }
        d2vec _tau_raw_i = orientation_class[i].get_tau_raw(distribution[indices[i][range(0, N)][range()]], shear_rate);
        _aggr_part = distribution[i][N][0];
        if (M != 1 && N != 1) { // We have orientation included
			for (int j = 0; j < 3; j++) {
				for (int k = 0; k < 3; k++) {
					_tau_raw[j][k] += _Npi * _aggr_part * _tau_raw_i[j][k];
				}
			}
		} else { //Only aggregation
			_tau_raw[s1_][s2_] += _Npi * _aggr_part;
		}
    }
    return _tau_raw[s1_][s2_];
}

///Get probability for each size class
void Combo_class::get_tot_prob(double time_now_) {
	double _aggr_prob = 0.0;

    cout << "ClassProb c = " << classes << " t = " << time_now_ << " ";

	for (int i = 0; i < classes; i++) {
        multi_array_ref<double, 2>::array_view<2>::type _tempView_wdist = distribution[indices[i][range(0, N)][range()]];
        cout << orientation_class[i].get_tot_prob(_tempView_wdist) << " ";
		_aggr_prob += distribution[i][N][0];
    }
    cout << _aggr_prob << " " << endl;
}

void Combo_class::save_aggr_distribution(double time_now_) {
    ofstream _save_file("aggr_dists.dat", std::ios::out | std::ios::app);
	_save_file << "#at " << time_now_ << endl;
	for (int i = 0; i < classes; i++) {
		_save_file << distribution[i][N][0] << endl;
	}
	_save_file << "\n\n";
}

void Combo_class::update_shear_rate(double shear_rate_) {
    shear_rate = shear_rate_;
    if (classes > 1) aggregation_class.update_v_aggregation(shear_rate_);
    if (M != 1 && N != 1) { //If this is false forget about the orientation part
        for (indexb i = 0; i < classes; i++) {
            orientation_class[i].initialize_abcde(shear_rate_);
            orientation_class[i].update_shear_rate(shear_rate_);
        }
    }
}