#include <iostream>
#include <cstdio>
#include <cmath>
#include <fstream>
#include <array>
using namespace std;

// Class to describe properties of a single discretized point in space (in cylindrical coord's).
// Contains spatial coordinates (x, r), viscosity (visc), velocity in x-direction (vx).
class point {
public:
    float x, r, visc, vx;

    // Function for printing values of a point, mostly for debugging.
    void print_values() const {
        printf("x=%f, r=%f, visc=%f, vx=%f\n", x, r, visc, vx);
    };
};

// Function to calculate fluid x-velocity in a single point in the pipe.
float vx_pipe(float dp, float visc, float l, float radius, float r) {
    return (dp/(4*visc*l))*(pow(radius, 2) - pow(r, 2));
}

// Writes r values and x-velocities along the radius to text file.
void write_r_vx(int nx, int ny, point radius[]) {
    ofstream file;
    file.open("../pyplot/r_vx.dat");
    for (int i = 0; i <= ny; i++) {
        file << radius[i].r << "\t" << radius[i].vx << "\n";
    }
    file.close();
}

int main() {
    float r, dp, visc, l, R;
    dp = 14.0; visc = 0.1; l = 0.2; R = 0.1;
    int nx, ny;
    ny = 5;

    point radius[ny];

    // Sets values for the point in the y-midpoint
    radius[0].r = 0;
    radius[0].visc = visc;
    radius[0].x = 1;
    radius[0].vx = vx_pipe(dp, visc, l, R, 0);
    // Loop goes through the other points in y-direction, starting from the middle
    for (int i = 1; i <= ny; i++) {
        r = ((float)i/ny)*R;
        radius[i].r = r;
        radius[i].visc = visc;
        radius[i].x = 1;
        radius[i].vx = vx_pipe(dp, visc, l, R, r);
    }

    for (int i = 0; i <= ny; i++) {
        radius[i].print_values();
    }
    write_r_vx(nx, ny, radius);

    return 0;
}
