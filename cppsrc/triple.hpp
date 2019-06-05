#ifndef TRIPLE_AO_2012
#define TRIPLE_AO_2012
#include <math.h>
#include <fstream>
#include <iostream> // remove this dependency from final version 
#ifndef PI
#define PI  3.14159265358979323846
#endif

///3d vector with conversion to spherical coordinates
///Notice -- you have to take care that the conversion is carried out
///by calling the appropiate transform function! 
class triple{
public:
    double x, y, z;
    double phi, theta, r;
    triple():x(0.0), y(0.0), z(0.0){;}
    triple(const double &phi_, const double &theta_){transform(phi_,theta_);}
    triple(const double &first_,const double &second_, const double &third_):x(first_),y(second_),z(third_){;}
    ///tranform() -- get "phi,theta,r" from current "x,y,z"
    triple& transform(){
         //~ r = sqrt(x*x + y*y + z*z);
        r = sqrt(pow(x,2) + pow(y,2) + pow(z,2));
        //~ if(r<1e-8){phi=0;theta=0;  
        if(r<1e-30){phi=0;theta=0;  
        //~ std::cerr << "x " << x << " y " << y << " z "  << z << std::endl;
        }
        else{
            phi = atan2(y,x);
            if(phi<0) phi+=2*PI;
            theta = acos(z/r);    
        }
        return *this;
    }
    ///transform(phi,theta) -- get "x,y,z" from "phi,theta" and r=1
    void transform(const double &phi_, const double &theta_){
        phi=phi_;
        theta=theta_;
        x=sin(theta_)*cos(phi_);
        y=sin(theta_)*sin(phi_);
        z=cos(theta_);
    }
    ///add_r(r) if you need also r setted -OBS! Assumes r now == 1. 
    triple& add_r(const double &r_){
        r=r_;
        x*=r;
        y*=r;
        z*=r;
        return *this;
    }
    ///transform_to_xyz() -- get "x,y,z" from current "phi,theta,r"
    triple& transform_to_xyz(){
        transform(phi,theta);
        add_r(r);
        return *this;
    }
    void shoutbox(int print_endl=1){
        std::cout <<x << " " << y << " " << z << " ";
        if(print_endl)
            std::cout << std::endl;
    }
    std::ofstream & shoutbox_file(std::ofstream &outfile ){
        outfile <<x << " " << y << " " << z << " "<<  std::endl;
        return outfile;
    }
    
    ///+= operator works in "xyz" but calls transform() afterwords 
    triple& operator +=(const triple& add_this){
        x+=add_this.x;
        y+=add_this.y;
        z+=add_this.z;
        transform();
        return *this;
    }
};
#endif
