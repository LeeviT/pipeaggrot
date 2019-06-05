#ifndef ORIENTATION_AO_2012
#define ORIENTATION_AO_2012
///The orientation class -- 3d orientation model (fokker-planck)
///Direct approach:: disrectization in spherical coordinates: control volume scheme
///Many functions defined in header since template functions ...

#include <vector>
#include "boost/multi_array.hpp"
#include "triple.hpp"
#include <fstream>
#include <iostream>
#ifndef PI
#define PI  3.14159265358979323846
#endif

typedef boost::multi_array<double, 2>  d2vec;
typedef boost::detail::multi_array::multi_array_view<double, 2> d2vec_view; 

class Orientation_ao{
private:
    int s1, s2;
    int runalready;
    ///Paramters of the system
    double  k_aspect, diffusion_coeff, shear_rate;
    ///The matrices 
    d2vec W_ij, W_ipj, W_imj, W_ijp, W_ijm, area;
    int N,M; ///N -- phi size and  M theta size
    std::vector<double> phi, theta;
    double h_phi, h_theta, tot_area;
    ///Set k_aspect from aspect ratio
    void set_k_aspect(const double aspect_ratio);
    ///Get rotation velocities according to shear rate (Jeffery equation)
    ///Currently this is only used in initialize_abcde, 
    ///so speed is not relevant
    triple dot_p(double phi_, double theta_);
    triple dot_p2(double phi_, double theta_);
    ///dot_p in spherical coordinates 
    std::pair<double, double> dot_p_angles(double phi_, double theta_);
    std::pair<double, double> dot_p_angles2(double phi_, double theta_);
    std::pair<double, double> dot_p_angles3(double phi_, double theta_);
    /// Initialize the matrices
    /// has to be called each time the shear rate changes (currently constant)
    /// PARAMETRIZE IF THIS HAS TO BE CALLED MORE OFTEN !!
    void initialize_abcde();
public:
    ///The main Initializer function
    void initialize(const double aspect_ratio, 
		double set_shear_rate, double set_diffusion_coeff,int M_,int N_, int s1_, int s2_);
    ///Return an approximately uniform distribution
    d2vec get_uniform_dist();
    ///Get dw_dt
    ///Template function here because of "array_view reference" != "array reference"
    ///Should now also eat vectors and arrays, because size information is not asked for
    template <class d2vec_view1>
    void get_dw_dt(d2vec_view1 &dw, const d2vec_view1 &w_old){ 
        int N_1=N-1, M_1=M-1;
        int i, j, im1, ip1, jm1,jp1;
        for(i=0;i!=N;++i){
            im1=i-1; ip1=i+1; if(!i) im1=N-1; else if(i==N_1) ip1=0;
            for(j=1; j!=M_1; j++){
                jm1=j-1; jp1=j+1;
                dw[i][j] += w_old[i][j] * W_ij[i][j] 
                    + w_old[i][jp1] * W_ijp[i][j]
                    + w_old[i][jm1] * W_ijm[i][j]
                    + w_old[ip1][j] * W_ipj[i][j] 
                    + w_old[im1][j] * W_imj[i][j];
            }
            //FOR control volume scheme:: flux over the poles is zero...
            //SPECIAL CASE:: NORTH POLE
            j=0; 
            dw[i][j] += w_old[i][j] * W_ij[i][j] 
                + w_old[i][1] * W_ijp[i][j] 
                + w_old[ip1][j] * W_ipj[i][j] 
                + w_old[im1][j] * W_imj[i][j];
            //SPECIAL CASE:: SOUTH POLE
            j=M_1;
            dw[i][j] += w_old[i][j]* W_ij[i][j] 
                + w_old[i][j-1] * W_ijm[i][j]
                + w_old[ip1][j] * W_ipj[i][j] 
                + w_old[im1][j] * W_imj[i][j];
        }
    }
    ///GET a2
    ///since we need w_dist, this is currently a template function ..
    template <class d2vec_view1>
    double get_a2_(const d2vec_view1 &w_dist, int index_1, int index_2){
        //get "a2" second order orientation tensor (function of w_dist)
        double a2[3][3]={0.0};
        for(int phi_now=0; phi_now!=N; phi_now++){
            for(int theta_now=0; theta_now!=M; theta_now++){
                triple p_triple(phi[phi_now],theta[theta_now]);  //These should be stored somewhere 
                double p_[3]={p_triple.x,p_triple.y,p_triple.z}; //or even these
                for(int i=0; i!=3; ++i){
                    for(int j=0; j!=3; ++j){
                        //basically one could even store something like a4_raw[phi_now][theta_now][i][j][k][l]
                        a2[i][j]+=p_[i]*p_[j]*area[phi_now][theta_now]*w_dist[phi_now][theta_now];
                    }
                }
            }
        }
        return a2[index_1][index_2];
    }
    
    ///A4 comp
    template <class d2vec_view1>
    double get_a4_comp(const d2vec_view1 &w_dist){
        //get "a2" second order orientation tensor (function of w_dist)
        double a2[3][3]={0.0};
        for(int phi_now=0; phi_now!=N; phi_now++){
            for(int theta_now=0; theta_now!=M; theta_now++){
                triple p_triple(phi[phi_now],theta[theta_now]);  //These should be stored somewhere 
                double p_[3]={p_triple.x,p_triple.y,p_triple.z}; //or even these
                for(int i=0; i!=3; ++i){
                    for(int j=0; j!=3; ++j){
                        //basically one could even store something like a4_raw[phi_now][theta_now][i][j][k][l]
                        a2[i][j]+=p_[i]*p_[j]*area[phi_now][theta_now]*w_dist[phi_now][theta_now];
                    }
                }
            }
        }
        double a4_comp [3][3][3][3]={0.0};
        double a4[3][3][3][3]={0.0};
        for(int phi_now=0; phi_now!=N; phi_now++){
            for(int theta_now=0; theta_now!=M; theta_now++){
                triple p_triple(phi[phi_now],theta[theta_now]);  //These should be stored somewhere 
                double p_[3]={p_triple.x,p_triple.y,p_triple.z}; //or even these
                for(int i=0; i!=3; ++i){
                    for(int j=0; j!=3; ++j){
                        for(int k=0; k!=3; ++k){
                            for(int l=0; l!=3; ++l){
                                //basically one could even store something like a4_raw[phi_now][theta_now][i][j][k][l]
                                a4[i][j][k][l]+=p_[i]*p_[j]*p_[k]*p_[l]*area[phi_now][theta_now]*w_dist[phi_now][theta_now];
                            }
                        }
                    }
                }
            }
        }
        //~ std::cout << "a4_0";
        for(int i=0; i!=3; ++i){
            for(int j=0; j!=3; ++j){
                for(int k=0; k!=3; ++k){
                    for(int l=0; l!=3; ++l){
                        a4_comp[i][j][k][l]=a2[i][j]*a2[k][l];
                        //~ std::cout << a4_comp[i][j][k][l] << " ";
                    }
                }
            }
        }
        //~ std::cout << sstd::endl;
        std::cout << "a4 \t a4_comp" << "\n";
        for(int i=0; i!=3; ++i){
            for(int j=0; j!=3; ++j){
                for(int k=0; k!=3; ++k){
                    for(int l=0; l!=3; ++l){
                        std::cout << a4[i][j][k][l] << "\t" << a4_comp[i][j][k][l] << 
                        "\t" << a4[i][j][k][l] - a4_comp[i][j][k][l] << "\n";
                    }
                }
            }
        }
        std::cout << std::endl;
        exit(1);
    }
    
    
    
    ///Return a4:D 
    ///TODO:: Optimize this heavily
    ///since we need w_dist, this is currently a template function ..
    template <class d2vec_view1>
    d2vec get_tau_raw(const d2vec_view1 &w_dist){
        d2vec a4_D(boost::extents[3][3]);
        std::fill(a4_D.data(),a4_D.data()+a4_D.num_elements(),0.0);
        //define D
        //D should be defined somewhere else
        double D[3][3]={0.0}; 
        D[s1][s2]=D[s2][s1]=shear_rate/2.0;
        //get "a4" fourth order orientation tensor (function of w_dist)
        double a4[3][3][3][3]={0.0};
        for(int phi_now=0; phi_now!=N; phi_now++){
            for(int theta_now=0; theta_now!=M; theta_now++){
                triple p_triple(phi[phi_now],theta[theta_now]);  //These should be stored somewhere 
                double p_[3]={p_triple.x,p_triple.y,p_triple.z}; //or even these                
                for(int i=0; i!=3; ++i){
                    for(int j=0; j!=3; ++j){
                        for(int k=0; k!=3; ++k){
                            for(int l=0; l!=3; ++l){
                                //basically one could even store something like a4_raw[phi_now][theta_now][i][j][k][l]
                                a4[i][j][k][l]+=p_[i]*p_[j]*p_[k]*p_[l]*area[phi_now][theta_now]*w_dist[phi_now][theta_now];
                            }
                        }
                    }
                }
            }
        }
        //Get a4:D
        for(int i=0; i!=3; ++i){
            for(int j=0; j!=3; ++j){
                for(int k=0; k!=3; ++k){
                    for(int l=0; l!=3; ++l){
                        a4_D[i][j]+=a4[i][j][k][l]*D[k][l];
                    }
                }
            }
        }
        return a4_D;
    }
    
    
    template <class d2vec_view1>
    ///Expression is:
    ///\sigma=-pI+2\mu(D+N_p(a_4:D))
    d2vec get_stress1(const d2vec_view1 &w_dist){
        std::cout << "GET STRESS " << std::endl;
        double Np=300; //TODO:: put this somewhere else
        //D should be defined somewhere else
        double D[3][3]={0.0}; 
        D[s1][s2]=D[s2][s1]=shear_rate/2.0;
        d2vec a4_D=get_tau_raw(w_dist);
        d2vec stress(boost::extents[3][3]);
        //~ double pressure=0;
        for(int i=0; i!=3;++i){
            for(int j=0; j!=3; ++j){
                stress[i][j]=2*(D[i][j]+Np*a4_D[i][j]);
                std::cout <<"#" << i << j << " " << stress[i][j]<< std::endl;

            }
            //~ stress[i][i]-=pressure;
        }
        return stress;
    }
    template <class d2vec_view1>
    double get_visc_part(d2vec_view1 &w_dist){
        normalize(w_dist); //w_dist is not const because of this line...
        double viscosity_part=0.0;
        double a4[3][3][3][3]={0.0};
        for(int phi_now=0; phi_now!=N; phi_now++){
            for(int theta_now=0; theta_now!=M; theta_now++){
                triple p_triple(phi[phi_now],theta[theta_now]);  //These should be stored somewhere 
                double p_[3]={p_triple.x,p_triple.y,p_triple.z}; //or even these
                for(int i=0; i!=3; ++i){
                    for(int j=0; j!=3; ++j){
                        for(int k=0; k!=3; ++k){
                            for(int l=0; l!=3; ++l){
                                a4[i][j][k][l]+=p_[i]*p_[j]*p_[k]*p_[l]*area[phi_now][theta_now]*w_dist[phi_now][theta_now];
                            }
                        }
                    }
                }
            }
        }
        //~ double a_xyxy = a4[0][1][0][1];
        double a_xyxy = a4[s1][s2][s1][s2];
        viscosity_part = 2*a_xyxy;
        return viscosity_part;
    }
    
   
    template <class d2vec_view1>
    void normalize(d2vec_view1 &w_dist){
    double tot_prob=0.0;
        for(int phi_now=0; phi_now!=N; phi_now++){
            for(int theta_now=0; theta_now!=M; theta_now++){
                 tot_prob += w_dist[phi_now][theta_now] * area[phi_now][theta_now];
            }
        }
        std::cout << "#tot_prob1= " << tot_prob << std::endl;
        //~ for(int phi_now=0; phi_now!=N; phi_now++){
            //~ for(int theta_now=0; theta_now!=M; theta_now++){
                //~ w_dist[phi_now][theta_now] /= tot_prob;
            //~ }
        //~ }
    }
    
    template <class d2vec_view1>
    double get_tot_prob(d2vec_view1 &w_dist){
        double tot_prob=0.0;
        for(int phi_now=0; phi_now!=N; phi_now++){
            for(int theta_now=0; theta_now!=M; theta_now++){
                 tot_prob += w_dist[phi_now][theta_now] * area[phi_now][theta_now];
            }
        }
        return tot_prob;
    }
    
    
    
    ///Get the "amount of particles" in this distribution 
    template <class d2vec_view1>
    double get_C(d2vec_view1 &w_dist){
        double tot_prob=0.0;
        for(int phi_now=0; phi_now!=N; phi_now++){
            for(int theta_now=0; theta_now!=M; theta_now++){
                 tot_prob += w_dist[phi_now][theta_now] * area[phi_now][theta_now];
            }
        }
        return tot_prob;
    }
    
    
    void area_test(){
        double tot_area=0.0;
        for(int phi_now=0; phi_now!=N; phi_now++){
            for(int theta_now=0; theta_now!=M; theta_now++){
             //~ std::cerr << theta_now << " " << phi_now << " " << area[phi_now][theta_now] << std::endl;
            tot_area += area[phi_now][theta_now];
            }
        }
        //~ std::cerr << tot_area/(4*PI) << std::endl;
    }
    
    ///Print the current distribution 
    ///TODO take the output stream as input parameter
    template <class d2vec_view1>
    //~ void shoutbox(const d2vec_view1 &w_dist, ostream & outfile=cout){
    void shoutbox(const d2vec_view1 &w_dist, int file_index){
        std::stringstream filename;
        filename << "class_" << file_index << ".dat";
        std::ofstream outfile;
        if(!runalready) outfile.open((filename.str()).c_str());
        else outfile.open((filename.str()).c_str(), std::ios::out | std::ios::app);
        runalready++;
        for(int i=0; i!=N; ++i){
            for(int j=0; j!=M; ++j){
                triple tmp(phi[i],theta[j]);
                //~ tmp.add_r(max(0.0,w_dist[i][j])/(4*PI));
                tmp.add_r(w_dist[i][j]/(4*PI));
                tmp.shoutbox_file(outfile);
            }
            outfile <<"\n";
        }
        int i=0;
        //Print first again in order to get grid in gnuplot
        for(int j=0; j!=M; ++j){
            triple tmp(phi[i],theta[j]);
            tmp.add_r(w_dist[i][j]/(4*PI));
            tmp.shoutbox_file(outfile);
        }
        outfile <<"\n";
        outfile << std::endl;
    }
    template <class d2vec_view1>
    double get_Np(const d2vec_view1 &dist_now){
        double Np_now=0.0;
        for(int i=0; i!=N; ++i){
            for(int j=0; j!=M; ++j){
               Np_now+=dist_now[i][j]*area[i][j];
            }
        }
        return Np_now;
    }
    ///Selfcheck -- mainly for Debugging
    void selfcheck();    
};
#endif
