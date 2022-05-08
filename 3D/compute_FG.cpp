#include <cmath>
using namespace std;

//math_diff.cpp
double du_dx(double ***var_u, int i, int j, int k, double dx);
double dv_dy(double ***var_v, int i, int j, int k, double dy);
double dw_dy(double ***var_w, int i, int j, int k, double dz);
double du2_dx(double ***var_u, int i, int j,int k, double dx, double gamma);
double dv2_dy(double ***var_v, int i, int j, int k, double dy, double gamma);
double dw2_dz(double ***var_w, int i, int j,int k, double dz, double gamma);
double duv_dx(double ***var_u, double ***var_v, int i, int j, int k, double dx, double gamma);
double duw_dz(double ***var_u, double ***var_w, int i, int j, int k, double dz, double gamma);
double duv_dy(double ***var_u, double ***var_v, int i, int j, int k, double dy, double gamma);
double dvw_dz(double ***var_v, double ***var_w, int i, int j, int k, double dz, double gamma);
double duw_dx(double ***var_u, double ***var_w, int i, int j, int k, double dx, double gamma);
double dvw_dy(double ***var_v, double ***var_w, int i, int j, int k, double dy, double gamma);
double d2u_dx2(double ***var_u, int i, int j, int k, double dx);
double d2v_dx2(double ***var_v, int i, int j, int k, double dx);
double d2w_dx2(double ***var_w, int i, int j, int k, double dx);
double d2u_dy2(double ***var_u, int i, int j, int k, double dy);
double d2v_dy2(double ***var_v, int i, int j, int k, double dy);
double d2w_dy2(double ***var_w, int i, int j, int k, double dy);
double d2u_dz2(double ***var_u, int i, int j, int k, double dz);
double d2v_dz2(double ***var_v, int i, int j, int k, double dz);
double d2w_dz2(double ***var_w, int i, int j, int k, double dz);


void compute_F(double ***var_F, double ***var_u, double ***var_v, double ***var_w, int nx, int ny, int nz, double dx, double dy, double dz, double dt, double gamma, double Re, double g_x) {
  for (int j = 1; j <= ny; j++) {
    for (int k = 1; k<= nz; k++){
      var_F[0][j][k] = var_u[0][j][k];
      var_F[nx][j][k] = var_u[nx][j][k];

      for (int i = 1; i <= nx - 1; i++) {
	var_F[i][j][k] = var_u[i][j][k] + dt*((d2u_dx2(var_u, i, j, k, dx) + d2u_dy2(var_u, i, j, k, dy)+ d2u_dz2(var_u, i, j, k, dz))/Re - du2_dx(var_u, i, j, k, dx, gamma) - duv_dy(var_u, var_v, i, j, k, dy, gamma) - duw_dz(var_u,var_w, i, j, k, dz, gamma) + g_x);
      }      
    }
  }
}
void compute_G(double ***var_G, double ***var_u, double ***var_v, double ***var_w, int nx, int ny, int nz, double dx, double dy, double dz, double dt, double gamma, double Re, double g_y) {
  for (int i = 1; i <= nx; i++) {
    for (int k = 1; k<=nz; k++){
      var_G[i][0][k] = var_v[i][0][k];
      var_G[i][ny][k] = var_v[i][ny][k];
      for (int j = 1; j <= ny-1; j++) {
	var_G[i][j][k] = var_v[i][j][k] + dt*( (d2v_dx2(var_v, i, j, k, dx) + d2v_dy2(var_v, i, j, k, dy) + d2v_dz2(var_v, i, j, k, dz))/Re - duv_dx(var_u, var_v, i, j, k, dx, gamma) - dv2_dy(var_v, i, j, k, dy, gamma) - dvw_dz(var_v, var_w, i, j, k, dz, gamma) + g_y);
      }
    }
  }
}

void compute_H(double ***var_H, double ***var_u, double ***var_v, double ***var_w, int nx, int ny, int nz, double dx, double dy, double dz, double dt, double gamma, double Re, double g_z) {
  for (int i = 1; i <= nx; i++) {
    for (int j = 1; j<=ny; j++){
      var_H[i][j][0] = var_w[i][j][0];
      var_H[i][j][nz] = var_w[i][j][nz];
      for (int k = 1; k <= nz-1; k++) {
	var_H[i][j][k] = var_w[i][j][k] + dt*( (d2w_dx2(var_w, i, j, k, dx) + d2w_dy2(var_w, i, j, k, dy) + d2w_dz2(var_w, i, j, k, dz))/Re - duw_dx(var_u, var_w, i, j, k, dx, gamma) - dw2_dz(var_w, i, j, k, dz, gamma) - dvw_dy(var_v, var_w, i, j, k, dy, gamma) + g_z);
      }
    }
  }
}


void compute_RHS(double ***var_RHS,double ***var_F, double ***var_G, double ***var_H, int nx, int ny, int nz, double dx, double dy, double dz, double dt) {
  for (int j = 1; j <= ny; j++) {
    for (int i = 1; i <= nx; i++) {
      for (int k = 1; k <= nz; k++){
	var_RHS[i][j][k] = ( (var_F[i][j][k]-var_F[i-1][j][k])/dx + (var_G[i][j][k]-var_G[i][j-1][k])/dy - (var_H[i][j][k]-var_H[i][j][k-1])/dz)/dt;
      }
    }
  }
}