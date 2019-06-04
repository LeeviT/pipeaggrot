#include <iostream>
#include <cstdio>
#include <cmath>
#include <fstream>
#include <array>
using namespace std;

// Class to describe properties of a single discretized point in space (in cylindrical coord's)
// Contains spatial coordinates (x, r), viscosity (visc), velocity in x-direction (vx)
class point {
public:
    float x, r, visc, vx;

    // Function for printing values of a point, mostly for debugging
    void print_values() const {
        printf("x=%f, r=%f, visc=%f, vx=%f\n", x, r, visc, vx);
    };
};


// Function to calculate fluid x-velocity in a single point in the pipe.
float vx_pipe(point single_point, float dp, float l, float radius) {
    return (dp/(4*single_point.visc*l))*(pow(radius, 2) - pow(single_point.r, 2));
}

// Writes r values and x-velocities along the radius to text file
void write_r_vx(int nx, int ny, point radius[]) {
    ofstream file;
    file.open("../pyplot/r_vx.dat");
    for (int i = 0; i <= ny; i++) {
        file << radius[i].r << "\t" << radius[i].vx << "\n";
    }
    file.close();
}

// Structure for handling input and returning it
struct input {
    float dp, l, visc, R;
    int nx, ny;
};

// Reads values from the input file to input structure and returns it
input read_input_file() {
    string a, b;
    float c;
    input values;
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
        } if (a == "nx") {
            values.nx = c;
        } if (a == "ny") {
            values.ny = c;
        }
    }
    file.close();
    return values;
}

int main() {
    // Variables for point class objects
    float r;
    // Initial values and constants for the system from the input file
    // dp=pressure difference, l=length of the pipe, R=radius of the pipe
    // nx=number of discretization points in x-direction, ny=same but in y-direction
    input input_params;
    float dp, l, visc, R;
    int nx, ny;

    // Assigns input values to local variables from input structure
    // Some values won't change during program runs so maybe could use values straight from struct, more error prone?
    input_params = read_input_file();
    dp = input_params.dp;
    l = input_params.l;
    visc = input_params.visc;
    R = input_params.R;
    nx = input_params.nx;
    ny = input_params.ny;
    printf("dp=%f, l=%f, visc=%f, R=%f, nx=%i, ny=%i\n", dp, l, visc, R, nx, ny);

    // Create point object, ny points along radius of the pipe
    point radius[ny];

    // Sets values for the point in the y-midpoint
    radius[0].r = 0;
    radius[0].visc = visc;
    radius[0].x = 1;
    radius[0].vx = vx_pipe(radius[0], dp, l, R);

    // Loop goes through the other points in y-direction, starting from the middle
    for (int i = 1; i <= ny; i++) {
        r = ((float)i/ny)*R;
        radius[i].r = r;
        radius[i].visc = visc;
        radius[i].x = 1;
        radius[i].vx = vx_pipe(radius[i], dp, l, R);
    }


    for (int i = 0; i <= ny; i++) {
        radius[i].print_values();
    }
    write_r_vx(nx, ny, radius);

    return 0;
}
