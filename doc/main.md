Documentation for the functions in **main.cpp** source code file.

* **int wfunc(realtype t_, N_Vector y_, N_Vector ydot_, void * userdata_):**  
  Function used by the SUNDIALS/CVODE solver in order to calculate weights or the ODE right-hand side
  for the solution.  
  _Arguments_: **t_** is the state variable, i.e. physical time in a current timestep, 
  **y_** is the current value of the dependent variable vector, i.e. current solution vector, **ydot_**
  is the output vector for _f(t, y)_, **userdata_** is the ```user_data``` pointer. Consult section
  **CVRhsFn** in CVODE documentation for more information.
  
* **vector<double> init_visc_vector(int ny_, double input_visc_):**
  Function to allocate memory for the viscosity vector and assign homogenous viscosity values to 
  vector.  
  _Arguments_: **ny_** number of discretization points along the radius, **input_visc_** input value for
  the viscosity, assumed to be homogenous everywhere in the beginning.
  
* **int main(int argc_, char * argv_[]):**
  The main function of the program. First initializes most variables and (class) instances used in the
  main function. Then initializes MPI, reads input file, prints input values and calculates how many
  discretization points each MPI rank handle. Later on, several vectors, e.g. viscosity vector, is 
  initialized, as well as Bessel function related stuff. After that, CVODE instances and options are 
  set for each discretization point. Once CVODE instances are ready, the main timestepping loop is
  started. Inside the timestepping loop, there's another loop going through discretization points. For
  each discretization radial coordinate r, viscosity, x coordinate are set and local flow velocity is
  calculated. Once local flow velocities/shear rates are calculated, they are updated to AO model. 
  However, this is not done for the point on the edge of the pipe since its velocity must be fixed to
  zero. Now CVODE time integrator is used advance solution in time and updated viscosity is obtained 
  from the updated solution. After obtaining new solution, program opens output files and writes
  viscosities, r values and flow velocities to output files. Once all the timesteps are done, program
  deallocates memory used by CVODE and destroys MPI instance.
  _Arguments_: **argc_** is a dummy argument for MPI, **argv_[]** is a command line argument vector for
  MPI.
