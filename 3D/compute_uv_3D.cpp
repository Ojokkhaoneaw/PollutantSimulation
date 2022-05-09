#include <cmath>
using namespace std;


void compute_uv_3D(double ***var_u, double ***var_v, double ***var_w, double ***var_F, double ***var_G, double ***var_H, double ***var_p_new, int nx, int ny, int nz, double dx, double dy, double dz, double dt) {

  for (int k = 1; k <= nz; k++) {
    for (int j = 1; j <= ny; j++) {
      for (int i = 1; i <= nx; i++) {
        if (i!= nx){
          var_u[i][j][k] = var_F[i][j][k] - (dt/dx)*(var_p_new[i+1][j][k]-var_p_new[i][j][k]);
        }
        if (j!=ny) {
          var_v[i][j][k] = var_G[i][j][k] - (dt/dy)*(var_p_new[i][j+1][k]-var_p_new[i][j][k]);
        }
        if (k!=nz) {
          var_w[i][j][k] = var_H[i][j][k] - (dt/dz)*(var_p_new[i][j][k+1]-var_p_new[i][j][k]);
	    }
      }
    }
  }
}
