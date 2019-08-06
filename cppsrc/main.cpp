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
#include "Combo_class.hpp"

#include "cvode/cvode.h"
#include "cvode/cvode_dense.h"
#include "cvode/cvode_lapack.h"
#include "nvector/nvector_serial.h"
#include "sundials/sundials_dense.h"
#include "sundials/sundials_types.h"
#include "sundials/sundials_math.h"
#include "sundials/sundials_spgmr.h"
#include "cvodes/cvodes_spgmr.h"
using namespace std;

// The function needed by the sundials solver "return dt"
// (get ydot at given t, y, and *userdata)
int wfunc(realtype t_, N_Vector y_, N_Vector ydot_, void * userdata_) {
    auto * _combined_ptr = (Combo_class *) userdata_;
    _combined_ptr->get_dt(NV_DATA_S(ydot_), NV_DATA_S(y_));
    return 0;
}

// Initialize viscosities containing vector which points radially from the midpoint outwards
// Viscosity is homogenous in the beginning, coming from the input file
vector<double> init_visc_vector(int ny_, double input_visc_) {
    vector<double> _visc;

    for (int i = 0; i <= ny_; i++) {
        _visc.push_back(input_visc_);
    }
    return _visc;
}

int main() {
    // Dummy variable for point class objects
    double _r;
    // Initialize variables and class instances
    input _params{};
    vector<Combo_class> _combined;
    vector<Combo_class *> _user_data;
    vector<double> _visc_vector;
    vector<point> _radius;
    double _time = 0, _shear_rate = 0.0, _upd_shear_rate, _visc_raw_til;
    int _system_size;
    // Initialize SUNDIALS stuff
    vector<void *> _cvode_mem;
    vector<N_Vector> _y0;
    vector<ofstream> _visctotfile;
    // Startup velocity stuff
    int _n_bessel_zeros = 200;

    // Read input file and assign input values to input::params struct
    _params.read_input_file();

    // Print input values to the screen
    printf("dp=%f, l=%f, visc=%f, R=%f, dt=%f, diff_coeff=%f\n", _params.getdp(), _params.getl(), _params.getvisc0(),
            _params.getR(), _params.getdt(), _params.getdiff_coeff());
    printf("aspect_ratio=%f, beta=%f, shear_rate_max=%f, Np_max=%f\n", _params.getaspect_ratio(), _params.getbeta(),
            _params.getshear_rate_max(), _params.getNp_max());
    printf("theta_grid=%i, phi_grid=%i, classes=%i, what_todo=%i\n", _params.gettheta_grid(), _params.getphi_grid(),
            _params.getclasses(), _params.getwhat_todo());
    printf("s1=%i, s2=%i, ny=%i, t0=%i, t_max=%i\n\n", _params.gets1(), _params.gets2(), _params.getny(),
            _params.gett0(), _params.gett_max());

    // Allocate memory and initialize uniform viscosity vector for the first timestep
    _visc_vector = init_visc_vector(_params.getny() + 1, _params.getvisc0());

    // Startup velocity stuff
    vector<double> _bessel_zeros(_n_bessel_zeros), _J_1_values(_n_bessel_zeros);
    vector<vector<double>> _J_0_values(_params.getny() + 1);
    for (int i = 0; i <= _params.getny(); i++) { _J_0_values[i].resize(_n_bessel_zeros); }
    _bessel_zeros = calc_bessel_zeros(_n_bessel_zeros);
    _J_0_values = J_0(_params.getny(), _bessel_zeros);
    _J_1_values = J_1(_bessel_zeros);

    // Initialize radius and combined vectors, ny+1 points along radius of the pipe
    // This must be done using *_back() function since vectors are accessed index-wise later on
    for (int i = 0; i <= _params.getny(); i++) {
        _radius.emplace_back();
        _combined.emplace_back();
        _visctotfile.emplace_back();
    }

    // Create ny AO-populations
    for (int i = 0; i <= _params.getny(); i++) {
        // Population system size in spherical coordinates
        _system_size = (_params.getphi_grid() + 1) * _params.gettheta_grid() * _params.getclasses();
        // Initialize AO model class instance (single population) using input values
        _combined.at(i).initialize(_shear_rate, _params);
        _user_data.push_back(&_combined.at(i));
        _y0.push_back(N_VNew_Serial(_system_size));
        // Set uniform distribution
        for (int j = 0; j < _system_size; j++) { _combined.at(i).set_uniform_state(NV_DATA_S(_y0.at(i))); }
        _cvode_mem.push_back(CVodeCreate(CV_ADAMS, CV_FUNCTIONAL));
        CVodeSetUserData(_cvode_mem.at(i), _user_data.at(i));
        // Initialize the integrator memory
        CVodeInit(_cvode_mem.at(i), wfunc, _params.gett0(), _y0.at(i));
        // Set the error weights
        CVodeSStolerances(_cvode_mem.at(i), 1.0e-8, 1.0e-8);
        CVodeSetMaxNumSteps(_cvode_mem.at(i), 100000);
        // Specify dense solver
        CVSpgmr(_cvode_mem.at(i), 0, 0);
        // We need to run the solver for very tiny timestep to initialize it properly:
        CVode(_cvode_mem.at(i), 1e-14, _y0.at(i), &_time, CV_ONE_STEP);
    }

    // The main timestep loop
    for (int t_step = 0; t_step <= _params.gett_max(); t_step++) {
        printf("t_step=%i of %i\n", t_step, _params.gett_max());
        cout << "time: " << _time << endl;

        // Loop goes through the other points in y-direction, starting from the middle and sets the following
        // values for each discretization point
        for (int i = 0; i <= _params.getny(); i++) {
            // If in y-midpoint of the pipe, r coordinate is zero, naturally
            _r = (i == 0) ? 0.0 : ((double) i / _params.getny()) * _params.getR();
            _radius[i].setr(_r);
            _radius[i].setvisc(_visc_vector[i]);
            _radius[i].setx(1.0);
            _radius[i].setvx(vx_pipe(_radius[i], _params, _bessel_zeros, _J_1_values, _J_0_values, i, t_step));

            // Update shear rate for the AO model, if not on the edge of the pipe
            if (i < _params.getny()) {
                _upd_shear_rate = (_radius[i].getvx() * 4.0) / _params.getR();
                _combined[i].update_shear_rate(_upd_shear_rate);
            }

            // Some AO model solver stuff
            if (i == 0) { CVode(_cvode_mem[i], _time + _params.getdt(), _y0[i], &_time, CV_NORMAL); }
            else { CVode(_cvode_mem[i], _time, _y0[i], &_time, CV_NORMAL); }
            _visc_raw_til = _combined[i].get_visc_raw(_params.gets1(), _params.gets2(), t_step);

            // Compute total viscosity of the fluid
            if (i < _params.getny()) {
                _visc_vector[i] = _params.getvisc0() + 2 * _params.getvisc0() * _visc_raw_til;
            } else {
                _visc_vector[i] = _params.getvisc0();
            }

            // Some printing stuff
            _visctotfile[i].open("../pysrc/visctot" + to_string(i) + ".dat", ios::app);
            _visctotfile[i] << _time << "\t" << _visc_vector[i] << endl;
            _visctotfile[i].close();
            cout << "r: " << _radius[i].getr() << ", total visc: " << _visc_vector[i] << ", v: "
                 << _radius[i].getvx() << endl;
            _combined[i].save_aggr_distribution(_time);
        }
        cout << " " << endl;
        // Writes r and vx values to file
        write_r_vx(_params.getny(), t_step, _radius);
    }

    for (int i = 0; i <= _params.getny(); i++) {
        _combined[i].save_state();
        N_VDestroy_Serial(_y0[i]);
        CVodeFree(&_cvode_mem[i]);
    }

    return 0;
}
