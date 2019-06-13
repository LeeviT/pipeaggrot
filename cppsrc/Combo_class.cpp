#include "Combo_class.hpp"
#include <algorithm>
#include <iostream>
#include <string.h> //memcpy is in string.h ...

typedef boost::multi_array_types::index_range range;
typedef d3vec::index indexb;
using namespace std;
using namespace boost;

///INITIALIZE
//Now we need extra space os size "classes" for the aggregation distribution
void Combo_class::initialize(double aspect_ratio, double shear_rate_, 
    double diffusion_coeff_, int classes_, int N_, int M_, 
    int s1, int s2, int selfcheck, double visco0_, double Np_max_, double beta_now, double shear_rate_max) {

    shear_rate=shear_rate_;
    diffusion_coeff_max=diffusion_coeff_; 
    classes=classes_; M=M_; N=N_;
    //Exra space of size classes*M for aggregation part <-- some overhead here...
    distribution.resize(extents[classes][N+1][M]);
    dist_dt.resize(extents[classes][N+1][M]);
    orientation_class.resize(classes); //This many orientation classes
    d_classes=1.0/classes;
    visco0=visco0_;
    Np_max=Np_max_;
    // By default aggregation model is not used but if condition is true, create aggregation class
    if(classes > 1)  aggregation_class.initialize(classes, shear_rate, beta_now, shear_rate_max);
    //Each size has its own orienation object
    //The diffusion coefficient is a function of the class
    if(M!=1 && N!=1){ //If this is false forget about the orientation part
        for(indexb i=0; i!=classes; ++i) {
	        double diffusion_coeff_i = diffusion_coeff_max*(i+1.0)/(classes*1.0);
            //~ cout << diffusion_coeff_i;  exit(1);
            orientation_class[i].initialize(aspect_ratio, shear_rate, diffusion_coeff_i, N, M, s1, s2);
        }
	}
    if(selfcheck==-1) {
        for(indexb i=0; i!=classes; ++i){
            orientation_class[i].selfcheck();
        }
        throw("selfcheck");
    }
}

///GET DT
///after this function the distribution stored in Combo_class is 
///the one given as input *y_now
void Combo_class::get_dt(double* y_dot, double* y_now)
{
    //Quite some overhead here since have to copy a large array at each step ... 
    //Working with boost multiarray -- that's why
    memcpy (distribution.data(), y_now, sizeof(double)*distribution.num_elements());
    //Initialize dist_dt with zeroes
    fill(dist_dt.data(),dist_dt.data()+dist_dt.num_elements(),0.0);
    //The tempView:s are created since need non-const reference to array view
    //Thus the array views cannot be created upon function call 
    //Probably there is an easier/cleaner way...
    ///Orientation part 
    if(M!=1 && N!=1) { //If this is flase forget about the orientation part
        for(int i=0; i!=classes; ++i){
            multi_array_ref<double, 2>::array_view<2>::type tempView_dt=dist_dt[indices[i][range(0,N)][range()]];
            orientation_class[i].get_dw_dt(tempView_dt, distribution[indices[i][range(0,N)][range()]]);
        }
	}
    ///Aggregation_fragmentation part
    if(classes>1)
    {
	    multi_array_ref<double, 1>::array_view<1>::type tempView_dt=dist_dt[indices[range()][N][0]];
        aggregation_class.get_wdist_dt(tempView_dt
		, distribution[indices[range()][N][0]]);
    }
    //Copy the calculated dist_dt (multiarray) to y_dot ...
    memcpy (y_dot, dist_dt.data(), sizeof(double)*dist_dt.num_elements());
}

///Simple function for setting uniform initial state
void Combo_class::set_uniform_state(double* y_now){
    //Step1:: the distribution
    fill(distribution.data(),distribution.data()+distribution.num_elements(),0.0);
    ifstream read_file("saved_state.dat");
    for (int i=0; i!=classes; ++i){
        for(int j=0; j!=N; ++j){
            for(int k=0; k!=M; ++k){
                distribution[i][j][k]= 1.0/(4.0*PI);
            }
        }
		distribution[i][N][0]= 1.0/classes;
    }
    //Step2:: copy the data from distribution to "the N_Vector"
    memcpy (y_now, distribution.data(), sizeof(double)*distribution.num_elements());
}

//Simple function for loading a state from file
void Combo_class::load_state(double* y_now){
    //Step1:: read the State from a file into distribution
    ifstream read_file("saved_state.dat");
    for (int i=0; i!=classes; ++i){
        for(int j=0; j!=N+1; ++j){
            for(int k=0; k!=M; ++k){
                read_file >> distribution[i][j][k];
            }
        }
    }
    //Step2:: copy the data from distribution to "the N_Vector"
    memcpy (y_now, distribution.data(), sizeof(double)*distribution.num_elements());
}
    
//Simple function for saving a state from a file
void Combo_class::save_state(){
    ofstream save_file("saved_state.dat");
    for (int i=0; i!=classes; ++i){
        for(int j=0; j!=N+1; ++j){
            for(int k=0; k!=M; ++k){
                save_file << distribution[i][j][k] << " ";
            }
        }
        save_file << "\n";
    }
}

double Combo_class::get_visc_raw(int s1_, int s2_) {
    double visc_raw=get_tau_raw(s1_,s2_);
    return (visc_raw/shear_rate);
}

///THIS should give the stress part due to the orientation...
double Combo_class::get_tau_raw(int s1_, int s2_){
    d2vec tau_raw(boost::extents[3][3]);
    std::fill(tau_raw.data(),tau_raw.data()+tau_raw.num_elements(),0.0); //fill with zeroes
    for(int i=0; i!=classes; ++i) {
        //Np varies linearly with the population variable n
        double Npi;
        if (classes==1) {
            Npi=Np_max;
        } else {
            Npi=Np_max*(i+1)*1.0/(classes);
        }
        d2vec tau_rawi = orientation_class[i].get_tau_raw(distribution[indices[i][range(0,N)][range()]]);
        double aggr_part=distribution[i][N][0];
        if(M!=1 && N!=1) { //We have orientation included
			for(int j=0; j!=3; ++j) {
				for(int k=0;k!=3; ++k) {
					tau_raw[j][k]+= Npi*aggr_part*tau_rawi[j][k];
				}
			}
		} else { //Only aggregation
			tau_raw[s1_][s2_]+= Npi*aggr_part;
		}
    }
    cout << "#TAURAW";
    for(int j=0; j!=3; ++j) {
        for(int k=0;k!=3; ++k) {
            cout << " " <<  tau_raw[j][k];
        }
    }
    cout << endl;
    //~ cout << tau_raw[s1_][s2_];
    return tau_raw[s1_][s2_];
}

///Get probability for each size class
void Combo_class::get_tot_prob(double time_now) {
    cout << "ClassProb c = " << classes << " t = " << time_now << " ";
	double aggr_prob = 0.0;
    for (int i = 0; i != classes; ++i) {
        multi_array_ref<double, 2>::array_view<2>::type tempView_wdist = distribution[indices[i][range(0,N)][range()]];
        cout << orientation_class[i].get_tot_prob(tempView_wdist) << " ";
		aggr_prob += distribution[i][N][0];
    }
    cout << aggr_prob << " " << endl;
}

void Combo_class::save_aggr_distribution(double time_now){
    ofstream save_file("aggr_dists.dat", std::ios::out | std::ios::app);
	save_file << "#at " << time_now << endl;
	for(int i=0; i!=classes; ++i){
		save_file << distribution[i][N][0] << endl;
	}
	save_file << "\n\n";
}

//
void Combo_class::update_shear_rate(double shear_rate_t) {
    shear_rate = shear_rate_t;
}