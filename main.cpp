#include <iostream>
#include <cmath>
using namespace std ;

void initialize(double **var, int nx, int ny) ;
void visualize(double **var, int nx, int ny) ;
void update(double **var, double **var_new, int nx, int ny) ;
void simulation(double **var, double **var_new, int nx, int ny, double c, double dx, double dy, double dt) ;




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
    // Boundary Condition
    phi[5][5] = 1. ;

    visualize(phi, nx, ny) ;
    for (int n = 0 ; n <= pow(10,1) ; n++){
        simulation(phi, phi_new,nx, ny, k, dx, dy, dt) ;
        update(phi, phi_new, nx, ny) ;
        visualize(phi, nx, ny) ;
    }
}
