#include <iostream>
#include <cstdio>
#include <cmath>
#include <fstream>
#include <array>
#include <ctime>
#include <cstdlib>
using namespace std;

// Class to describe properties of a single discretized point in space (in cylindrical coord's)
// Contains spatial coordinates (x, r), viscosity (visc), velocity in x-direction (vx)
class point {
private:
    float x, r, visc, vx;

public:
    // Function for printing values of a point, mostly for debugging
    void print_values() const {
        printf("x=%f, r=%f, visc=%f, vx=%f\n", x, r, visc, vx);
    };
    // Setter functions for each variable
    void setx(float x_set) { x = x_set; };
    void setr(float r_set) { r = r_set; };
    void setvisc(float visc_set) { visc = visc_set; };
    void setvx(float vx_set) { vx = vx_set; }
    // Getter functions for each variable
    float getx() { return x; };
    float getr() { return r; };
    float getvisc() { return visc; };
    float getvx() { return vx; };
};


// Function to calculate fluid x-velocity in a single point in the pipe
float vx_pipe(point single_point, float dp, float l, float radius) {
    return (dp/(4*single_point.getvisc()*l))*(pow(radius, 2) - pow(single_point.getr(), 2));
}

// Writes r values and x-velocities along the radius to text file
void write_r_vx(int ny, int t_step, point radius[]) {
    ofstream file;
    // string filename = "../pysrc/r_vx" + to_string(t_step) + ".dat";
    file.open("../pysrc/r_vx" + to_string(t_step) + ".dat");
    for (int i = 0; i <= ny; i++) {
        file << radius[i].getr() << "\t" << radius[i].getvx() << "\n";
    }
    file.close();
}

// Structure for handling input and returning it
struct input {
    float dp, l, visc, R;
    int ny, t_max;
};

// Reads values from the input file to input structure and returns it
input read_input_file() {
    string a, b;
    float c;
    input values{};
    ifstream file;
    file.open("input.txt");

    while (file >> a >> b >> c) {
        if (a == "pressure_diff") {
            values.dp = c;
        } if (a == "pipe_length") {
            values.l = c;
        } if (a == "visc") {
            values.visc = c;
        } if (a == "radius") {
            values.R = c;
        } if (a == "ny") {
            values.ny = (int)c;
        } if (a == "t_max") {
            values.t_max = (int)c;
        }
    }
    file.close();
    return values;
}

// Initialize viscosities containing vector which points radially from the midpoint outwards
float *init_visc_vector(int ny, float input_visc) {
    float *visc = new float[ny];
    float tmprand;

    srand (static_cast <unsigned> (time(nullptr)));

    for (int i = 0; i <= ny; i++) {
        tmprand = static_cast <float> (rand()) / (static_cast <float> (RAND_MAX/0.1)) + 1.0;
        visc[i] = input_visc*tmprand;
    }
    return visc;
}

// Sets modified viscosity values to viscosity vector each timestep
float *set_visc_vector(int ny, float *visc_vector) {
    float tmprand;

    srand (static_cast <unsigned> (time(nullptr)));

    for (int i = 0; i <= ny; i++) {
        tmprand = static_cast <float> (rand()) / (static_cast <float> (RAND_MAX/0.02)) + 1.0;
        visc_vector[i] *= tmprand;
    }
    return visc_vector;
}

int main() {
    // Variables for point class objects
    float r;
    // Initial values and constants for the system from the input file
    // dp=pressure difference, l=length of the pipe, R=radius of the pipe
    // nx=number of discretization points in x-direction, ny=same but in y-direction, t_max=max number of timesteps
    input input_params{};
    float dp, l, visc, R;
    float *visc_vector;
    int ny, t_max;

    // Assigns input values to local variables using input structure
    // Some values won't change during program runs so maybe could use values straight from struct, more error prone?
    input_params = read_input_file();
    dp = input_params.dp;
    l = input_params.l;
    visc = input_params.visc;
    R = input_params.R;
    ny = input_params.ny;
    t_max = input_params.t_max;
    printf("dp=%f, l=%f, visc=%f, R=%f, ny=%i, t_max=%i\n", dp, l, visc, R, ny, t_max);

    // Allocate memory and initialize viscosity vector for the first timestep
    visc_vector = init_visc_vector(ny, visc);

    // Create point object, ny points along radius of the pipe
    point radius[ny];

    for (int t_step = 0; t_step <= t_max; t_step++) {
        printf("t_step=%i\n", t_step);
        // Set new viscosity values if not in the very first timestep
        if (t_step > 0) {
            set_visc_vector(ny, visc_vector);
        }

        // Sets values for the point in the y-midpoint
        radius[0].setr(0.0);
        radius[0].setvisc(visc_vector[0]);
        radius[0].setx(1.0);
        radius[0].setvx(vx_pipe(radius[0], dp, l, R));

        // Loop goes through the other points in y-direction, starting from the middle
        for (int i = 1; i <= ny; i++) {
            r = ((float) i / ny) * R;
            radius[i].setr(r);
            radius[i].setvisc(visc_vector[i]);
            radius[i].setx(1.0);
            radius[i].setvx(vx_pipe(radius[i], dp, l, R));
        }

        // Prints values along the radius
        for (int i = 0; i <= ny; i++) {
            radius[i].print_values();
        }

        // Writes r and values to file
        write_r_vx(ny, t_step, radius);
    }

    delete [] visc_vector;

    return 0;
}
