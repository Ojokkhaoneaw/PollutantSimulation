#include <cmath>
using namespace std;

//math_diff.cpp
double du_dx(double **var_u, int i, int j, double dx);
double dv_dy(double **var_v, int i, int j, double dy);
double du2_dx(double **var_u, int i, int j, double dx, double gamma);
double dv2_dy(double **var_v, int i, int j, double dy, double gamma);
double duv_dx(double **var_u, double **var_v, int i, int j, double dx, double gamma);
double duv_dy(double **var_u, double **var_v, int i, int j, double dy, double gamma);
double d2u_dx2(double **var_u, int i, int j, double dx);
double d2v_dx2(double **var_v, int i, int j, double dx);
double d2u_dy2(double **var_u, int i, int j, double dy);
double d2v_dy2(double **var_v, int i, int j, double dy);

void compute_F(double **var_F, double **var_u, double **var_v, int nx, int ny, double dx, double dy, double dt, double gamma, double Re, double g_x) {
  for (int j = 1; j <= ny; j++) {
    
    var_F[0][j] = var_u[0][j];
    var_F[nx][j] = var_u[nx][j];
    
    for (int i = 1; i <= nx - 1; i++) {
      var_F[i][j] = var_u[i][j] + dt*((d2u_dx2(var_u, i, j, dx) + d2u_dy2(var_u, i, j, dy))/Re - du2_dx(var_u, i, j, dx, gamma) - duv_dy(var_u, var_v, i, j, dy, gamma) + g_x);
    }      
  }
}

void compute_G(double **var_G, double **var_u, double **var_v, int nx, int ny, double dx, double dy, double dt, double gamma, double Re, double g_y) {
  for (int i = 1; i <= nx; i++) {
    
    var_G[i][0] = var_v[i][0];
    var_G[i][ny] = var_v[i][ny];
    
    for (int j = 1; j <= ny-1; j++) {
      
      var_G[i][j] = var_v[i][j] + dt*( (d2v_dx2(var_v, i, j, dx) + d2v_dy2(var_v, i, j, dy))/Re - duv_dx(var_u, var_v, i, j, dx, gamma) - dv2_dy(var_v, i, j, dy, gamma) + g_y);

    }    
  }
}

void compute_RHS(double **var_RHS,double **var_F, double **var_G, int nx, int ny, double dx, double dy, double dt) {
  for (int j = 1; j <= ny; j++) {
    for (int i = 1; i <= nx; i++) {
      
      var_RHS[i][j] = ( (var_F[i][j]-var_F[i-1][j])/dx + (var_G[i][j]-var_G[i][j-1])/dy )/dt;

    }
  }
}
