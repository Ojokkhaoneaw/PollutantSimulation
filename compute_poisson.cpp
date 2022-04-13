#include <cmath> 
#include <iostream>
using namespace std;

void boundary_pressure(double **var_p, int nx, int ny, double dx, double dy, char wall, char type, double P);

void update(double **var, double **var_new, int nx, int ny);

void poisson(double **var_p, double **var_p_new, double **RHS, int nx, int ny, double dx, double dy, double omega, double eps, int iter_max) {

  int iter = 0;
  double n, e, s, w;
  double r_norm;
  
  while (r_norm<eps && iter<iter_max){
    
    double r_ij = 0;
    iter += 1;
    
    for (int i = 1; i <= nx; i++) {

      var_p_new[i][0] = var_p[i][1];
      var_p_new[i][ny+1] = var_p[i][ny];
	
      if (i==1) { e=1; w=0; }
      else if (i==nx) { e=0; w=1; }
      else { e = 1; w = 1; }

      for (int j = 1; j <= ny; j++) {
	
	var_p_new[0][j] = var_p[1][j];
	var_p_new[nx+1][j] = var_p[nx][j];
	
	if (j == 1) { n=1; s=0; }
	else if (j == ny) { n = 0; s = 1; }
	else { n=1; s=1; }
	var_p_new[i][j] = (1-omega)*var_p[i][j] + omega/( ((e+w)/pow(dx, 2.))+((n+s)/pow(dy, 2.)) )
	  * ( (e*var_p[i+1][j]+w*var_p_new[i-1][j])/pow(dx, 2.) + (n*var_p[i][j+1]+s*var_p_new[i][j-1])/pow(dy, 2.) - RHS[i][j] );

	r_ij += pow(( e*(var_p[i+1][j]-var_p[i][j])-w*(var_p[i][j]-var_p[i-1][j]) )/pow(dx, 2.) + ( n*(var_p[i][j+1]-var_p[i][j])-s*(var_p[i][j]-var_p[i][j-1]) )/pow(dy, 2.) - RHS[i][j], 2.);
	
      }
    }
    // Update & BC Presure
    update(var_p, var_p_new, nx, ny);
    boundary_pressure(var_p, nx, ny, dx, dy, 'north', 'N', 0);
    boundary_pressure(var_p, nx, ny, dx, dy, 'east', 'D', 0);
    boundary_pressure(var_p, nx, ny, dx, dy, 'south', 'N', 0);
    boundary_pressure(var_p, nx, ny, dx, dy, 'west', 'N', 0);
    
    r_norm = sqrt(r_ij/(nx*ny));

    
    }  

}
