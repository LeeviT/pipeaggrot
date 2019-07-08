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
//~ typedef boost::multi_array<double, 2>  d2vec;
//~ typedef boost::detail::multi_array::multi_array_view<double, 2> d2vec_view; 


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
void Orientation_ao::set_k_aspect(const double aspect_ratio) {
    double aspect_two = pow(aspect_ratio, 2);
    k_aspect = (aspect_two - 1.0)/(aspect_two + 1.0);
}

//dot_p in spherical coordinates 
//--do not return nans put -1e-5 instead -- careful here...
//first=dot_p_theta -- second=dot__phi
pair<double, double> Orientation_ao::dot_p_angles(double phi_, double theta_, double shear_rate_t) {
    //x y z 
    triple r_;
    r_.x = sin(theta_)*cos(phi_);
    r_.y = sin(theta_)*sin(phi_);
    r_.z = cos(theta_);
    //phi unit vector
    triple p_;
    p_.x = -sin(phi_);
    p_.y = cos(phi_);
    p_.z = 0;
    //theta unit vector
    triple t_;
    t_.x = cos(theta_)*cos(phi_);
    t_.y = cos(theta_)*sin(phi_);
    t_.z = -sin(theta_);
    pair<double,double> dot_p_s;
    double kappa[3][3] = {{0.0, 0.0, 0.0,}, {0.0,0.0,0.0}, {0.0,0.0,0.0}};
    double kappa_t[3][3] = {{0.0, 0.0, 0.0,}, {0.0,0.0,0.0}, {0.0,0.0,0.0}};
    kappa[s2][s1] = kappa_t[s1][s2] = shear_rate_t;
    double lambda = k_aspect;
    double lambda_p = (lambda + 1)/2.0;
    // now lambda m==0
    double lambda_m = (lambda - 1)/2.0;
    double k_rp = doublepoint_operator(kappa, r_, p_);
    double kt_rp = doublepoint_operator(kappa_t, r_, p_);
    double k_rt = doublepoint_operator(kappa, r_, t_);
    double kt_rt = doublepoint_operator(kappa_t, r_, t_);
    dot_p_s.second  = lambda_m*kt_rp + lambda_p*k_rp;
    dot_p_s.first  = lambda_m*kt_rt + lambda_p*k_rt;
    return dot_p_s;
}
//----------------------------------------------------------------------
//----------------------------------------------------------------------
//----------------------------------------------------------------------
///PUBLIC FUNCTIONS::
//initialize these -- 
// has to be called each time the shear rate changes (currently constant)
// PARAMETRIZE IF THIS HAS TO BE CALLED MORE OFTEN !!
// Since Only thing that actually changes 
// upon shear rate modification is  dot_p_angles
void Orientation_ao::initialize_abcde(double shear_rate_t_) {
    //in W_???: "p" refers to +1 and "m" to i-1  (imj->i-1,j) etc.
    W_ij.resize(extents[N][M]);
    W_ijp.resize(extents[N][M]);
    W_ijm.resize(extents[N][M]);
    W_ipj.resize(extents[N][M]);
    W_imj.resize(extents[N][M]);
    area.resize(extents[N][M]);
    tot_area = 0.0;
    double phi_, phi_m, phi_p, theta_, theta_m, theta_p;
    for (int i = 0; i < N; i++) {
        phi_= phi[i];
        phi_p = phi[i] + 0.5*h_phi;
        phi_m = phi[i] - 0.5*h_phi;
        for (int j = 0; j < M; j++) {
            theta_ = theta[j];
            theta_p = theta[j] + 0.5*h_theta;
            theta_m = theta[j] - 0.5*h_theta;
            //Area is exactly given by this formulation:
            area[i][j] = h_phi*fabs(cos(theta_m) - cos(theta_p));
            if (theta_m < PI/2 && theta_p > PI/2)
                area[i][j] = h_phi*(2-(1-fabs(cos(theta_m))) - (1-fabs(cos(theta_p))));
                //~ area[i][j]=h_phi*(fabs(cos(theta_m))+fabs(cos(theta_p)));
            W_ipj[i][j] =   
                (-sin(theta_)/2*dot_p_angles(phi_p, theta_, shear_rate_t_).second / sin(theta_) //orientation part
                + diffusion_coeff / (sin(theta_)*h_phi) )  //diffusion part
                //~ /(sin(theta_)*h_phi);
                * (h_theta/area[i][j]);
            W_imj[i][j] = 
                ( +sin(theta_)/ 2*dot_p_angles(phi_m, theta_, shear_rate_t_).second / sin(theta_) //orientation part
                + diffusion_coeff / (sin(theta_)*h_phi) )  //diffusion part
                //~ /(sin(theta_)*h_phi);
                *(h_theta/area[i][j]);
            W_ijp[i][j] =  
                ( -sin(theta_p)/ 2*dot_p_angles(phi_, theta_, shear_rate_t_).first //orientation part
                + diffusion_coeff * sin(theta_p) / h_theta ) //diffusion part
                //~ /(sin(theta_)*h_theta);
                *(h_phi/area[i][j]);
            W_ijm[i][j] = 
                ( +sin(theta_m)/2*dot_p_angles(phi_,theta_m, shear_rate_t_).first
                + diffusion_coeff * sin(theta_m) / h_theta )
                //~ /(sin(theta_)*h_theta);
                *(h_phi/area[i][j]);
            W_ij[i][j] = W_ijp[i][j] +  W_ijm[i][j] + W_ipj[i][j] + W_imj[i][j]
                - 2*diffusion_coeff * sin(theta_m) * h_phi / (h_theta * area[i][j]) 
                - 2*diffusion_coeff * sin(theta_p) * h_phi / (h_theta * area[i][j])
                - 4*diffusion_coeff * h_theta / (sin(theta_)*h_phi * area[i][j]);
                //~ - 2*diffusion_coeff * sin(theta_m) / (sin(theta_)*h_theta*h_theta)
                //~ - 2*diffusion_coeff * sin(theta_p) / (sin(theta_)*h_theta*h_theta)
                //~ - 4*diffusion_coeff / (sin(theta_)*sin(theta_)*h_phi*h_phi);
            //~ area[i][j] = h_theta*sin(theta_)*h_phi; 
            tot_area += area[i][j];
        }
    }
    //~ cout << "TOT area " << tot_area << endl; exit(1);
}

void Orientation_ao::initialize(const double aspect_ratio, double set_shear_rate, double set_diffusion_coeff,
                                int N_, int M_, int s1_, int s2_) {
    runalready = 0; //this is for deciding whether to create new files or append
    s1 = s1_;
    s2 = s2_;
    M = M_;
    N = N_; //system size
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
    shear_rate = set_shear_rate; //shear rate
    set_k_aspect(aspect_ratio); //aspect ratio initialization
    diffusion_coeff = set_diffusion_coeff; //diffusion coefficient
    initialize_abcde(set_shear_rate);
}
//Return an approximately uniform distribution
d2vec Orientation_ao::get_uniform_dist() {
    d2vec w_uniform(extents[N][M]);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            w_uniform[i][j] = 1.0/(tot_area);
        }
    }
    return w_uniform;
}
//Get dw_dt --- in header...
//Template function here because of "array_view reference" != "array reference"
//Functions defined in header ...


//----------------------------------------------------------------------
//----------------------------------------------------------------------
//----------------------------------------------------------------------
///Selfcheck -- mainly for Debugging
void Orientation_ao::selfcheck() {
    cout << "\n\n";
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            cout << (phi[i]/(2*PI))*360 << " " << (theta[j]/(2*PI))*360 << " "
                 << (dot_p_angles(phi[i], theta[j], shear_rate)).second << endl;
        }
    }
    cout << "\n\n";
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            cout << (phi[i]/(2*PI))*360 << " " << (theta[j]/(2*PI))*360 << " "
                 << (dot_p_angles(phi[i], theta[j], shear_rate)).first << endl;
        }
    }
    cout << "\n\n";
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            cout << (phi[i]/(2*PI))*360 << " " << (theta[j]/(2*PI))*360 << " " <<  W_ij[i][j] << endl;
        }
    }
    cout << "\n\n";
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            cout << (phi[i]/(2*PI))*360 << " " << (theta[j]/(2*PI))*360 << " " << W_ijp[i][j] << endl;
        }
    }
    cout << "\n\n";
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            cout << (phi[i]/(2*PI))*360 << " " << (theta[j]/(2*PI))*360 << " " << W_ijm[i][j]<< endl;
        }
    }
    cout << "\n\n";
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            cout << (phi[i]/(2*PI))*360 << " " << (theta[j]/(2*PI))*360 << " " << W_ipj[i][j] << endl;
        }
    }
    cout << "\n\n";
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            cout << (phi[i]/(2*PI))*360 << " " << (theta[j]/(2*PI))*360 << " " <<  W_imj[i][j] << endl;
        }
    }
    cout << "\n\n";
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            double dt_unity = W_ij[i][j] + W_ijp[i][j] + W_ijm[i][j] + W_ipj[i][j] + W_imj[i][j];
            cout << (phi[i]/(2*PI))*360 << " " << (theta[j]/(2*PI))*360 << " " << dt_unity << endl;
        }
    }
    cout << "\n\n";
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            triple temp(phi[i],theta[j]);
            double dt_unity = W_ij[i][j] + W_ijp[i][j] + W_ijm[i][j]; //+W_ipj[i][j]+ W_imj[i][j];
            cout << temp.x << " " << temp.y <<" " << temp.z << " " << dt_unity << endl;
        }
    }
    cout << "\n\n";

    double dt_unity = 0;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            dt_unity = W_ij[i][j] + W_ijp[i][j] + W_ijm[i][j] + W_ipj[i][j] + W_imj[i][j];
            if (j == 0) {
                triple temp1(phi[i], theta[j]);
                temp1.add_r(dt_unity);
                cout << temp1.x << "\t" << temp1.y << "\t" <<  temp1.z << endl;
            } else if (j == M) {
                triple temp2(phi[i], theta[j]);
                temp2.add_r(dt_unity);
                cout << temp2.x << "\t" << temp2.y << "\t" <<  temp2.z << endl;
            } else {
                triple temp(phi[i], theta[j]);
                temp.add_r(dt_unity);
                temp.transform();
                cout << temp.x << "\t" << temp.y << "\t" << temp.z << endl;
            }
        }
    }
    cout << "\n\n";
    
    //Check dot_p_angles -- 2 cases
    cout << "\n\n";
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            triple temp(phi[i], theta[j]);
            pair<double, double> temp_dt_angles = dot_p_angles(phi[i], theta[j], shear_rate);
            temp.add_r(temp_dt_angles.first);
            cout << temp.x << "\t" << temp.y << "\t"<<  temp.z << endl;
        }             
    }
    cout << "\n\n";
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            triple temp(phi[i], theta[j]);
            pair<double, double> temp_dt_angles = dot_p_angles(phi[i], theta[j], shear_rate);
            temp.add_r(temp_dt_angles.second);
            cout << temp.x << "\t" << temp.y << "\t"<<  temp.z << endl;
        }             
    }
    cout << "\n\n";
    //Check dot_p recreate -- 3 cases
    cout << "\n\n";
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            triple temp(phi[i], theta[j]);
            pair<double, double> temp_dt_angles = dot_p_angles(phi[i], theta[j], shear_rate);
            double px = cos(theta[j])*cos(phi[i])*temp_dt_angles.first
            - sin(theta[j])*sin(phi[i])*temp_dt_angles.second; 
            temp.add_r(px);
            cout << temp.x << "\t" << temp.y << "\t"<<  temp.z << endl;
        }             
    }
    cout << "\n\n";
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            triple temp(phi[i], theta[j]);
            pair<double, double> temp_dt_angles = dot_p_angles(phi[i], theta[j], shear_rate);
            double py = cos(theta[j])*sin(phi[i])*temp_dt_angles.first
            + sin(theta[j])*cos(phi[i])*temp_dt_angles.second; 
            temp.add_r(py);
            cout << temp.x << "\t" << temp.y << "\t"<<  temp.z << endl;
        }             
    }
    cout << "\n\n";
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            triple temp(phi[i], theta[j]);
            pair<double, double> temp_dt_angles = dot_p_angles(phi[i], theta[j], shear_rate);
            double pz = -sin(theta[j])*temp_dt_angles.first;
            temp.add_r(pz);
            cout << temp.x << "\t" << temp.y << "\t"<<  temp.z << endl;
        }             
    }
    cout << "\n\n";
}
