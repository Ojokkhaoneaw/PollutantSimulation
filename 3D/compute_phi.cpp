#include <cmath>
using namespace std;
double d2phi_dx2(double ***var_phi, int i, int j, int k, double dx);
double d2phi_dy2(double ***var_phi, int i, int j, int k, double dy);
double d2phi_dz2(double ***var_phi, int i, int j, int k, double dz);

double duphi_dx(double ***var_u, double ***var_phi, int i, int j, int k, double dx, double gamma) ;
double dvphi_dy(double ***var_v, double ***var_phi, int i, int j, int k, double dy, double gamma) ;
double dwphi_dz(double ***var_w, double ***var_phi, int i, int j, int k, double dx, double gamma) {


  void compute_phi(double ***var_phi,double ***var_phinew,double ***var_u, double ***var_v, double ***var_w, int nx, int ny, int nz, double dx, double dy, double dz, double dt,double gamma, double Re){

  for (int i = 1; i<=nx; i++){
    for (int j = 1; j<=ny; j++){
      for (int k = 1; k<=nz; k++){
	var_phinew[i][j][k] = var_phi[i][j][k] + dt*((d2phi_dx2(var_phi, i, j, k, dx) + d2phi_dy2(var_phi, i, j, k, dy) + d2phi_dz2(var_phi, i, j, k, dz))/Re - ( duphi_dx(var_u,var_phi, i, j, k,dx,gamma) + dvphi_dy(var_v,var_phi,i ,j , k,dy,gamma) + dwphi_dz(var_w,var_phi,i ,j , k,dz,gamma) )) ;
      }
    }
  }


}
