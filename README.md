# pipeaggrot

## General description
Program combining non-Newtonian pipe flow and aggregation-orientation (AO) model. A fluid flowing 
in a pipe is thixotropic -- in this case nanocellulose suspension. Particles in suspension are 
rigid and ellipsoidal. In order to find out rheological properties of suspension, behavior of 
particles must be calculated first. To do this, program is using particle distributions instead of 
single particles in order to calculate distribution function for particle aggregations and 
orientations. From AO distribution, local stress of fluid is calculated. Basically, in order to get
stress, program is calculating 
![](doc/stress.png)  
where the fourth moment of orientation tensor is calculated from
![](doc/moment.png)  
See e.g. paper [Rheological modeling of carbon nanotube aggregate suspensions 
(W. K. A. Ma et al. 2008)
](https://sor.scitation.org/doi/abs/10.1122/1.2982932) for more detailed description of AO model. 
As a local stress of fluid is calculated, it is used to calculate a local viscosity from which 
fluid velocity is calculated in discretization points in a following way
![](doc/velo.png)  
where λ<sub>n</sub> are roots of Bessel function of the first kind of the first order, i.e. 
J<sub>0</sub>(λ<sub>n</sub>)=0.  

## Source code
Detailed description of source code functions (those written by me and few others) can be found 
under [```pipeaggrot/doc/```](https://github.com/LeeviT/pipeaggrot/tree/master_par/doc) directory.

## Installing 
### Requirements 
Program requires C++11 and any compiler which supports it. Currently in the makefile Intel C++ 
compiler is used.  
Also [Boost C++ Libraries](https://www.boost.org/) and [SUNDIALS/CVODE ODE solver
library](https://computing.llnl.gov/projects/sundials/cvode) are required. Installation 
instructions of those libraries can be found in their documentation. However, both libraries must 
be installed under [```pipeaggrot/cppsrc/```](https://github.com/LeeviT/pipeaggrot/tree/master_par/cppsrc) 
directory. So, if using ```cmake``` in SUNDIALS installation, use the following commands  
```
cmake -DCMAKE_INSTALL_PREFIX=$SOME_PATH/pipeaggrot/cppsrc \
-DEXAMPLES_INSTALL_PATH=$SOME_PATH/pipeaggrot/cppsrc/examples \
-DOPENMP_ENABLE=ON \
../sundials-VERSION
```
