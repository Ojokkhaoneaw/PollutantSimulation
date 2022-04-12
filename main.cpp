#include <iostream>
#include <cmath>
using namespace std ;

//void initialize(double **var, int nx, int ny) ;
//void visualize(double **var, int nx, int ny) ;
//void update(double **var, double **var_new, int nx, int ny) ;
//void simulation(double **var, double **var_new, int nx, int ny, double c, double dx, double dy, double dt) ;
// Initialize variable
void initialize(double **var, int nx, int ny) {
    for (int i = 0 ; i <= nx-1 ; i++){
        for (int j = 0 ; j <= ny-1 ; j++) {
            var[i][j] = 0.0 ;
        }
    }
}

// Visualize variable
void visualize(double **var, int nx, int ny) {
    for (int i = 0 ; i <= nx-1 ; i++){
        for (int j = 0 ; j <= ny-1 ; j++) {
            cout << var[i][j] << " " ;
        }
        cout << "\n" ; ;
    }
}

//Update variable
void update(double **var, double **var_new, int nx, int ny) {
    for (int i = 0 ; i <= nx+1 ; i++){
        for (int j = 0 ; j <= ny+1 ; j++) {
            var[i][j] = var_new[i][j] ;
        }
    }
}

// Simulate 2D heat equation
void simulation(double **var, double **var_new, int nx, int ny, double c, double dx, double dy, double dt){
    double RHS ;
    // Simulation
    for (int i = 1 ; i <= nx-2 ; i++) {
        for(int j = 1 ; j <= ny-2; j++) {
            RHS = c * ((var[i +1][j] - 2.*var[i][j] + var[i-1][j]) / (dx* dy))+
                  c * ((var[i][j+1] - 2.*var[i][j] + var[i][j-1]) / (dx* dy)) ;
            var_new[i][j] = var[i][j] + dt * RHS;
        }
    }
}

int main() {
    int nx = 11 ;
    int ny = 11 ;
    double k = 1. ;
    double dx = 1. ;
    double dy = 1. ;
    double dt = 0.01 ;
//    --------------------------------
    double **phi ;
    phi = (double **) malloc (nx * sizeof(double)) ;
    for (int i = 0 ; i < nx ; i++) {
        phi[i] = (double *) malloc (ny * sizeof(double))  ;
    }
    double **phi_new ;
    phi_new = (double **) malloc (nx * sizeof(double)) ;
    for (int i = 0 ; i < nx ; i++) {
        phi_new[i] = (double *) malloc (ny * sizeof(double))  ;
    }
//    --------------------------------
    initialize(phi, nx, ny) ;
    initialize(phi_new, nx, ny) ;
    cout << "HELLO" << "\n" ;
    for (int n = 0 ; n <= pow(10,1) ; n++) {
        simulation(phi, phi_new,nx, ny, k, dx, dy, dt) ;
        update(phi, phi_new, nx, ny) ;
        visualize(phi, nx, ny) ;
        cout << n << "\n" ;
    }
}
