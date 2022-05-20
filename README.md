# Pollutant Simulation
 Final project of Applied Numerical Computation in Engineering
### 2D Simulation :
1. run `g++ -o run main.cpp math_diff.cpp compute_uv.cpp compute_poisson.cpp compute_uv.cpp compute_phi.cpp && ./run`
2. `main.cpp` is a main file for 2D simulation 
3. The other files are function files. 
   + `Libraries.cpp`  contains 
     + Boundary condition function
     + Initialize-variable function
     + Update-variable function
     + Save-file function 
     + Read-file function
     + Visualize-variable function 
     + Export-vtk-file function
     + Export-dat-file function
   + `math_diff.cpp` is used for computing derivative variables.
   + `compute_FG.cpp` is used for computing F and G term in momentum equations.
   + `compute_poisson.cpp` is used for computing pressure at the next time step
   + `compute_uv.cpp` is used for computing x-velocity and y-velocity at the next time step
   + `compute_phi.cpp` is used for computing pollutant at the next time step.  


### 3D Simulation
1. run `g++ -o run main_3D.cpp math_diff_3D.cpp compute_uv_3D.cpp compute_poisson_3D.cpp compute_uv_3D.cpp compute_phi_3D.cpp && ./run`
2. `main_3D.cpp` is a main file for 2D simulation.
3. The other files are function files.
    + `Libraries_3D.cpp`  contains
        + Boundary condition function
        + Initialize-variable function
        + Update-variable function
        + Save-file function
        + Read-file function
        + Visualize-variable function
        + Export-vtk-file function
        + Export-dat-file function
    + `math_diff_3D.cpp` is used for computing derivative variables.
    + `compute_FG_3D.cpp` is used for computing F and G term in momentum equations.
    + `compute_poisson_3D.cpp` is used for computing pressure at the next time step
    + `compute_uv_3D.cpp` is used for computing x-velocity, y-velocity and z-velocity at the next time step
    + `compute_phi_3D.cpp` is used for computing pollutant at the next time step. 

### Results
+  `Link` : https://drive.google.com/drive/folders/1jikQWdbfAfHTeM6nLEkycZDPEzcpUCtx?usp=sharing