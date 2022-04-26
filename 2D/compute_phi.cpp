#include <cmath>
using namespace std;
double d2phi_dx2(double **var_phi, int i, int j, double dx);
double d2phi_dy2(double **var_phi, int i, int j, double dy);
double duphi_dx(double **var_u, double **var_phi, int i, int j, double dx, double gamma) ;
double dvphi_dy(double **var_v, double **var_phi, int i, int j, double dy, double gamma) ;

void compute_phi(double **var_phi,double **var_phinew,double **var_u, double **var_v,  int nx, int ny,double dx,double dy,double dt,double gamma, double Re){

  for (int i = 1; i<=nx; i++){
    for (int j = 1; j<=ny; j++){

      var_phinew[i][j] = var_phi[i][j] + dt*((d2phi_dx2(var_phi, i, j, dx) + d2phi_dy2(var_phi, i, j, dy))/Re - ( duphi_dx(var_u,var_phi,i ,j ,dx,gamma) + dvphi_dy(var_v,var_phi,i ,j ,dy,gamma) )) ;

    }
  }


}
