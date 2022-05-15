#include <cmath> 
#include <iostream>
using namespace std;

void pressure_condition_3D(double ***var_p, int nx, int ny, int nz, double dx, double dy, double dz, char wall, char type, double P);

void update_3D(double ***var, double ***var_new, int nx, int ny, int nz);

void poisson_3D(double ***var_p, double ***var_p_new, double ***RHS, int nx, int ny, int nz, double dx, double dy, double dz, double omega, double eps, int iter_max) {

  int iter = 0;
  double n, e, s, w, b, f;
  double r_norm = 1.;
  
  while (r_norm > eps && iter<iter_max){
    
    double r_ijk = 0;
    iter += 1;
    
    for (int i = 1; i <= nx; i++) {
        for (int j = 1; j <= ny; j++) {
	        for (int k = 1; k <= nz; k++){
               var_p_new[i][0][k] = var_p[i][1][k];
               var_p_new[i][ny+1][k] = var_p[i][ny][k];

               var_p_new[0][j][k] = var_p[1][j][k];
               var_p_new[nx+1][j][k] = var_p[nx][j][k];

               var_p_new[i][j][0] = var_p[1][j][1];
               var_p_new[i][j][nz+1] = var_p[nx][j][nz];

               if (i==1) { e=1; w=0; }
               else if (i==nx) { e=0; w=1; }
               else { e = 1; w = 1; }

               if (j == 1) { n=1; s=0; }
               else if (j == ny) { n = 0; s = 1; }
               else { n=1; s=1; }

               if (k == 1) { f=1; b=0; }
               else if (k == nz) { f = 0; b = 1; }
               else { f=1; b=1; }

               var_p_new[i][j][k] = (1-omega)*var_p[i][j][k] + omega/( ((e+w)/pow(dx, 2.))+((n+s)/pow(dy, 2.))+((b+f)/pow(dz, 2.)) ) *
               ((e*var_p[i+1][j][k]+w*var_p_new[i-1][j][k])/pow(dx, 2.) +
                (n*var_p[i][j+1][k]+s*var_p_new[i][j-1][k])/pow(dy, 2.) +
                (f*var_p[i][j][k+1]+b*var_p_new[i][j][k-1])/pow(dz, 2.) -
                RHS[i][j][k] );

	            r_ijk += pow(( e*(var_p[i+1][j][k]-var_p[i][j][k])-w*(var_p[i][j][k]-var_p[i-1][j][k]) )/pow(dx, 2.) + ( n*(var_p[i][j+1][k]-var_p[i][j][k])-s*(var_p[i][j][k]-var_p[i][j-1][k]) )/pow(dy, 2.) + ( f*(var_p[i][j][k+1]-var_p[i][j][k])-b*(var_p[i][j][k]-var_p[i][j][k-1]) )/pow(dz, 2.) - RHS[i][j][k], 2.);
	        }
       }
    }
    // Update & BC Presure
    update_3D(var_p, var_p_new, nx, ny, nz);
    pressure_condition_3D(var_p, nx, ny, nz, dx, dy, dz, 'n', 'N', 0);
    pressure_condition_3D(var_p, nx, ny, nz, dx, dy, dz, 'e', 'N', 0);
    pressure_condition_3D(var_p, nx, ny, nz, dx, dy, dz, 's', 'N', 0);
    pressure_condition_3D(var_p, nx, ny, nz, dx, dy, dz, 'w', 'N', 0);
    pressure_condition_3D(var_p, nx, ny, nz, dx, dy, dz, 'f', 'N', 0);
    pressure_condition_3D(var_p, nx, ny, nz, dx, dy, dz, 'b', 'N', 0);
    r_norm = sqrt(r_ijk/(nx*ny*nz));
  }  
}
