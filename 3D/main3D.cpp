#include <iostream>
#include <cmath>
using namespace std ;

void initialize3D(double ***var, int nx, int ny ,int nz, double c) ;
void visualize(double **var, int nx, int ny) ;
void visualize3D(double ***var, int nx, int ny, int nz) ;
void update(double **var, double **var_new, int nx, int ny) ;
void paraview(int num_iter, const string& varName, double **var, int nx, int ny, double dx, double dy) ;
void paraview3D(int num_iter, const string& varName, double ***var, int nx, int ny, int nz, double dx, double dy,double dz) ;
void noslip_condition(double **var_u, double **var_v, int nx, int ny, char wall);
void outflow_condition(double **var_u, double **var_v, int nx, int ny, char wall);
void inflow_condition(double **var_u, double **var_v, int nx, int ny, char wall, double u, double v);
void pressure_condition(double **var_P, int nx, int ny, double dx, double dy, char wall, char type, double P) ;
void phi_condition(double **var_phi,int nx,int ny,double dx,double dy, char side, double phi) ;


void compute_F(double **var_F, double **var_u, double **var_v, int nx, int ny, double dx, double dy, double dt, double gamma, double Re, double g_x);
void compute_G(double **var_G, double **var_u, double **var_v, int nx, int ny, double dx, double dy, double dt, double gamma, double Re, double g_y);
void compute_RHS(double **var_RHS,double **var_F, double **var_G, int nx, int ny, double dx, double dy, double dt) ;


void poisson(double **var_p, double **var_p_new, double **RHS, int nx, int ny, double dx, double dy, double omega, double eps, int iter_max);


void compute_uv(double **var_u, double **var_v, double **var_F, double **var_G, double **var_p_new, int nx, int ny, double dx, double dy, double dt);
void compute_phi(double **var_phi,double **var_phinew,double **var_u, double **var_v,  int nx, int ny,double dx,double dy,double dt,double gamma, double Re) ;


int main() {
    const int iter(5000);
    const int iter_max(1000);
    const int nx(10); // 100
    const int ny(10); // 10
    const int nz(10) ;
    const double dy = 1;
    const double dx = 1.*dy ;
    const double dz = 1.*dy ;// 3
    const double dt = 0.01;
    const double u_init = 1.;
    const double v_init = 0.;
    const double p_init = 0.;
    const double phi_init = 1. ;
    const double Re = 150.; // Reynolds number
    const double g_x = 0.;
    const double g_y = 0.;
    const double g_z = 0. ;
    const double gamma = 0.5; // Upwind differencing factor
    const double omega = 1.7; // Relaxation parameter for SOR iteration
    const double eps = 0.00001; // Stopping tolerance for pressure iteration
    int i , j, k  ;

    // ------------------------------------------------
    double ***u; // u^(n)
    u = (double ***) malloc ((nx + 2) * sizeof(double));
    for (i = 0; i <= nx + 1; i++) {
        u[i] = (double **) malloc((ny + 2) * sizeof(double));
        for (j = 0; j <= ny + 1; j++) {
            u[i][j] = (double *) malloc((nz + 2) * sizeof(double));
        }
    }
    // ------------------------------------------------
    double ***v; // u^(n)
    v = (double ***) malloc ((nx + 2) * sizeof(double));
    for (i = 0; i <= nx + 1; i++) {
        v[i] = (double **) malloc((ny + 2) * sizeof(double));
        for (j = 0; j <= ny + 1; j++) {
            v[i][j] = (double *) malloc((nz + 2) * sizeof(double));
        }
    }
    // ------------------------------------------------
    double ***p; // u^(n)
    p = (double ***) malloc ((nx + 2) * sizeof(double));
    for (i = 0; i <= nx + 1; i++) {
        p[i] = (double **) malloc((ny + 2) * sizeof(double));
        for (j = 0; j <= ny + 1; j++) {
            p[i][j] = (double *) malloc((nz + 2) * sizeof(double));
        }
    }
    // ------------------------------------------------

    double ***F; // u^(n)
    F = (double ***) malloc ((nx + 2) * sizeof(double));
    for (i = 0; i <= nx + 1; i++) {
        F[i] = (double **) malloc((ny + 2) * sizeof(double));
        for (j = 0; j <= ny + 1; j++) {
            F[i][j] = (double *) malloc((nz + 2) * sizeof(double));
        }
    }
    // ------------------------------------------------

    double ***G; // u^(n)
    G = (double ***) malloc ((nx + 2) * sizeof(double));
    for (i = 0; i <= nx + 1; i++) {
        G[i] = (double **) malloc((ny + 2) * sizeof(double));
        for (j = 0; j <= ny + 1; j++) {
            G[i][j] = (double *) malloc((nz + 2) * sizeof(double));
        }
    }
    // ------------------------------------------------

    double ***RHS; // u^(n)
    RHS = (double ***) malloc ((nx + 2) * sizeof(double));
    for (i = 0; i <= nx + 1; i++) {
        RHS[i] = (double **) malloc((ny + 2) * sizeof(double));
        for (j = 0; j <= ny + 1; j++) {
            RHS[i][j] = (double *) malloc((nz + 2) * sizeof(double));
        }
    }
    // ------------------------------------------------

    double ***p_new; // u^(n)
    p_new = (double ***) malloc ((nx + 2) * sizeof(double));
    for (i = 0; i <= nx + 1; i++) {
        p_new[i] = (double **) malloc((ny + 2) * sizeof(double));
        for (j = 0; j <= ny + 1; j++) {
            p_new[i][j] = (double *) malloc((nz + 2) * sizeof(double));
        }
    }
    // ------------------------------------------------

    double ***phi; // u^(n)
    phi = (double ***) malloc ((nx + 2) * sizeof(double));
    for (i = 0; i <= nx + 1; i++) {
        phi[i] = (double **) malloc((ny + 2) * sizeof(double));
        for (j = 0; j <= ny + 1; j++) {
            phi[i][j] = (double *) malloc((nz + 2) * sizeof(double));
        }
    }
    // ------------------------------------------------

    double ***phi_new; // u^(n)
    phi_new = (double ***) malloc ((nx + 2) * sizeof(double));
    for (i = 0; i <= nx + 1; i++) {
        phi_new[i] = (double **) malloc((ny + 2) * sizeof(double));
        for (j = 0; j <= ny + 1; j++) {
            phi_new[i][j] = (double *) malloc((nz + 2) * sizeof(double));
        }
    }
    // ------------------------------------------------

    initialize3D(u,nx,ny,nz,1) ;
    initialize3D(v,nx,ny,nz,v_init) ;
    initialize3D(p,nx,ny,nz,p_init) ;
    initialize3D(phi,nx,ny,nz,0) ;
    for (i = 0 ; i <= nx+1 ; i++){
        for (j = 0 ; j <= ny+1 ; j++) {
            for (k = 0 ; k <= nz+1 ; k++) {
                u[i][j][k] = i*i+j*j+k*k ;
            }
        }
    }
    visualize3D(u,nx,ny,nz) ;
    paraview3D(1, "u", u, nx, ny, nz, dx, dy,dz) ;
//    // West
//    inflow_condition(u, v, nx, ny, 'w', u_init, v_init);
//    pressure_condition(p, nx, ny, dx, dy, 'w', 'N', p_init);
//    phi_condition(phi,nx,ny,dx,dy,'w',phi_init) ;
//    // East
//    outflow_condition(u, v, nx, ny, 'e');
//    pressure_condition(p, nx, ny, dx, dy, 'e', 'D', p_init);
//    phi_condition(phi,nx,ny,dx,dy,'e',phi_init) ;
//    // Wall (Top)
//    noslip_condition(u, v, nx, ny, 'n');
//    pressure_condition(p, nx, ny, dx, dy, 'n', 'N', p_init);
//    phi_condition(phi,nx,ny,dx,dy,'n',phi_init) ;
//    // Wall (South)
//    noslip_condition(u, v, nx, ny, 's');
//    pressure_condition(p, nx, ny, dx, dy, 's', 'N', p_init);
//    phi_condition(phi,nx,ny,dx,dy,'s',phi_init) ;
//    for (int num_iter = 1; num_iter <= iter; num_iter++) {
//
//        compute_F(F, u, v, nx, ny, dx, dy, dt, gamma, Re, g_x);
//        compute_G(G, u, v, nx, ny, dx, dy, dt, gamma,Re, g_y);
//        compute_RHS(RHS, F, G, nx, ny, dx, dy, dt);
//
//        poisson(p, p_new, RHS, nx, ny, dx, dy, omega, eps, iter_max);
//        inflow_condition(u ,v , nx, ny, 'w', u_init, v_init);
//        outflow_condition(u ,v , nx, ny, 'e');
//        noslip_condition(u ,v , nx, ny, 'n');
//        noslip_condition(u ,v , nx, ny, 's');
//
//        compute_uv(u, v, F, G, p_new, nx, ny, dx, dy, dt);
//        compute_phi(phi,phi_new,u,v,nx,ny,dx,dy,dt,gamma,Re) ;
//        update(phi,phi_new,nx,ny) ;
//        phi_condition(phi,nx,ny,dx,dy,'w',phi_init) ;
//        phi_condition(phi,nx,ny,dx,dy,'e',phi_init) ;
//        phi_condition(phi,nx,ny,dx,dy,'n',phi_init) ;
//        phi_condition(phi,nx,ny,dx,dy,'s',phi_init) ;
//        if (num_iter == 1 || num_iter%5 == 0) {
//            cout << "time step : " << num_iter <<"\n" ;
//            cout << "------------phi-------------" <<"\n" ;
//            visualize(phi,nx,ny);
//            paraview(num_iter ,"p" , p, nx, ny, dx, dy);
//            paraview(num_iter ,"F" , F, nx, ny, dx, dy);
//            paraview(num_iter ,"u" , u, nx, ny, dx, dy);
//            paraview(num_iter ,"v" , v, nx, ny, dx, dy) ;
//            paraview(num_iter ,"phi" , phi, nx, ny, dx, dy);
//        }

//    }
}