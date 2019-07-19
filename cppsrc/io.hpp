#ifndef CPPSRC_IO_HPP
#define CPPSRC_IO_HPP

using namespace std;

// Structure for handling input and returning it
struct input {
private:
    // beta = v_d/v_c
    // dp=pressure difference, l=length of the pipe, R=radius of the pipe
    // nx=number of discretization points in x-direction, ny=same but in y-direction, t_max=max number of timesteps
    double dp, l, visc0, R, dt, diff_coeff, aspect_ratio, beta, shear_rate_max, Np_max;
    int theta_grid, phi_grid, classes, what_todo, s1, s2, ny, t0, t_max;
public:
    // Setter functions for each variable
    void setdp(double dp_);
    void setl(double l_);
    void setvisc0(double visc0_);
    void setR(double R_);
    void setdt(double dt_);
    void setdiff_coeff(double diff_coeff_);
    void setaspect_ratio(double aspect_ratio_);
    void setbeta(double beta_);
    void setshear_rate_max(double shear_rate_max_);
    void setNp_max(double Np_max_);
    void settheta_grid(int theta_grid_);
    void setphi_grid(int phi_grid_);
    void setclasses(int classes_);
    void setwhat_todo(int what_todo_);
    void sets1(int s1_);
    void sets2(int s2_);
    void setny(int ny_);
    void sett0(int t0_);
    void sett_max(int t_max_);

    // Getter functions for each variable
    double getdp() const;
    double getl() const;
    double getvisc0() const;
    double getR() const;
    double getdt() const;
    double getdiff_coeff() const;
    double getaspect_ratio() const;
    double getbeta() const;
    double getshear_rate_max() const;
    double getNp_max() const;
    int gettheta_grid() const;
    int getphi_grid() const;
    int getclasses() const;
    int getwhat_todo() const;
    int gets1() const;
    int gets2() const;
    int getny() const;
    int gett0() const;
    int gett_max() const;

    // Function reading input values from the input file
    void read_input_file();
};

// Function writing r and x-velo values to file
void write_r_vx(int ny_, int t_step_, vector<point> radius_);

#endif