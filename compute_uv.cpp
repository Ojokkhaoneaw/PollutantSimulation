#include <cmath>
using namespace std;


void compute_uv(double **var_u, double **var_v, double **var_F, double **var_G, double **var_p_new, int nx, int ny, double dx, double dy, double dt) {

  
  for (int j = 1; j <= ny; j++) {
    for (int i = 1; i <= nx; i++) {
      if (i!= nx){
	    var_u[i][j] = var_F[i][j] - (dt/dx)*(var_p_new[i+1][j]-var_p_new[i][j]);
      }
      if (j!=ny) {
	    var_v[i][j] = var_G[i][j] - (dt/dy)*(var_p_new[i][j+1]-var_p_new[i][j]);
      }
     
    }
  } 
  
}
