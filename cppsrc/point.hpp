#ifndef CPPSRC_POINT_HPP
#define CPPSRC_POINT_HPP

// Class to describe properties of a single discretized point in space (in cylindrical coord's)
// Contains spatial coordinates (x, r), viscosity (visc), velocity in x-direction (vx)
class point {
private:
    double x, r, visc, vx;
public:
    // Constructor
    point();

    // Function for printing values of a point, mostly for debugging
    void print_values();

    // Setter functions for each variable
    void setx(double x);
    void setr(double r);
    void setvisc(double visc);
    void setvx(double vx);

    // Getter functions for each variable
    double getx();
    double getr();
    double getvisc();
    double getvx();
};

#endif
