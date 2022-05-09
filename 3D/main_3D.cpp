#include <iostream>
#include <cmath>
using namespace std ;

void initialize_3D(double ***var, int nx, int ny ,int nz, double c) ;
//void visualize_3D(double ***var, int nx, int ny, int nz) ;
void update_3D(double ***var, double ***var_new, int nx, int ny, int nz) ;
void paraview_3D(int num_iter, const string& varName, double ***var, int nx, int ny, int nz, double dx, double dy,double dz) ;
void noslip_condition_3D(double ***var_u, double ***var_v,double ***var_w, int nx, int ny, int nz, char side);
void outflow_condition_3D(double ***var_u, double ***var_v,double ***var_w, int nx, int ny, int nz , char side) ;
void inflow_condition_3D(double ***var_u, double ***var_v,double ***var_w, int nx, int ny, int nz, char side, double u, double v, double w);
void pressure_condition_3D(double ***var_P, int nx, int ny, int nz,  double dx, double dy, double dz, char side, char type, double P) ;
//void phi_condition_3D(double ***var_phi,int nx,int ny,double dx,double dy, char side, double phi) ;

void compute_F_3D(double ***var_F, double ***var_u, double ***var_v, double ***var_w, int nx, int ny, int nz, double dx, double dy, double dz, double dt, double gamma, double Re, double g_x);
void compute_G_3D(double ***var_G, double ***var_u, double ***var_v, double ***var_w, int nx, int ny, int nz, double dx, double dy, double dz, double dt, double gamma, double Re, double g_y);
void compute_H_3D(double ***var_H, double ***var_u, double ***var_v, double ***var_w, int nx, int ny, int nz, double dx, double dy, double dz, double dt, double gamma, double Re, double g_z);
void compute_RHS_3D(double ***var_RHS,double ***var_F, double ***var_G, double ***var_H, int nx, int ny, int nz, double dx, double dy, double dz, double dt);

void poisson_3D(double ***var_p, double ***var_p_new, double ***RHS, int nx, int ny, int nz, double dx, double dy, double dz, double omega, double eps, int iter_max);


void compute_uv_3D(double ***var_u, double ***var_v, double ***var_w, double ***var_F, double ***var_G, double ***var_H, double ***var_p_new, int nx, int ny, int nz, double dx, double dy, double dz, double dt) ;
//void compute_phi_3D(double ***var_phi,double ***var_phinew,double ***var_u, double ***var_v, double ***var_w, int nx, int ny, int nz, double dx, double dy, double dz, double dt,double gamma, double Re) ;


int main() {
    const int iter(100);
    const int iter_max(1000);
    const int nx(100); // 100
    const int ny(9); // 10
    const int nz(9) ;
    const double dy = 0.1;
    const double dx = 2.*dy ;
    const double dz = 2.*dy ;// 3
    const double dt = 0.01;
    const double u_init = 1.;
    const double v_init = 0.;
    const double w_init = 0.;
    const double p_init = 0.;
    const double phi_init = 1. ;
    const double Re = 150.; // Reynolds number
    const double g_x = 0.;
    const double g_y = 0.;
    const double g_z = 0.;
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
    double ***w; // u^(n)
    w = (double ***) malloc ((nx + 2) * sizeof(double));
    for (i = 0; i <= nx + 1; i++) {
        w[i] = (double **) malloc((ny + 2) * sizeof(double));
        for (j = 0; j <= ny + 1; j++) {
            w[i][j] = (double *) malloc((nz + 2) * sizeof(double));
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
    double ***H; // u^(n)
    H = (double ***) malloc ((nx + 2) * sizeof(double));
    for (i = 0; i <= nx + 1; i++) {
        H[i] = (double **) malloc((ny + 2) * sizeof(double));
        for (j = 0; j <= ny + 1; j++) {
            H[i][j] = (double *) malloc((nz + 2) * sizeof(double));
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

    initialize_3D(u,nx,ny,nz,1) ;
    initialize_3D(v,nx,ny,nz,v_init) ;
    initialize_3D(w,nx,ny,nz,v_init) ;
    initialize_3D(p,nx,ny,nz,p_init) ;
    // West
    inflow_condition_3D(u, v, w, nx, ny, nz, 'w', u_init, v_init, w_init);
    pressure_condition_3D(p, nx, ny, nz, dx, dy, dz, 'w', 'N', p_init);
//    phi_condition_3D(phi,nx,ny,dx,dy,'w',phi_init) ;
    // East
    outflow_condition_3D(u, v, w, nx, ny, nz, 'e');
    pressure_condition_3D(p, nx, ny, nz, dx, dy, dz, 'e', 'D', p_init);
//    phi_condition_3D(phi,nx,ny,dx,dy,'e',phi_init) ;
    // Wall (Top)
    noslip_condition_3D(u, v, w, nx, ny, nz, 'n');
    pressure_condition_3D(p, nx, ny, nz, dx, dy, dz, 'n', 'N', p_init);
//    phi_condition_3D(phi,nx,ny,dx,dy,'n',phi_init) ;
    // Wall (South)
    noslip_condition_3D(u, v,w, nx, ny, nz, 's');
    pressure_condition_3D(p, nx, ny, nz, dx, dy, dz, 's', 'N', p_init);
//    phi_condition_3D(phi,nx,ny,dx,dy,'s',phi_init) ;
    // Wall (Front)
    noslip_condition_3D(u, v,w, nx, ny, nz, 'f');
    pressure_condition_3D(p, nx, ny, nz, dx, dy, dz, 'f', 'N', p_init);
//    phi_condition_3D(phi,nx,ny,dx,dy,'f',phi_init) ;
    // Wall (Behind)
    noslip_condition_3D(u, v,w, nx, ny, nz, 'b');
    pressure_condition_3D(p, nx, ny, nz, dx, dy, dz, 'b', 'N', p_init);
//    phi_condition_3D(phi,nx,ny,dx,dy,'b',phi_init) ;
    for (int num_iter = 1; num_iter <= iter; num_iter++) {

        compute_F_3D(F, u, v, w, nx, ny, nz, dx, dy, dz, dt, gamma, Re, g_x);
        compute_G_3D(G, u, v, w, nx, ny, nz, dx, dy, dz, dt, gamma, Re, g_y);
        compute_H_3D(H, u, v, w, nx, ny, nz, dx, dy, dz, dt, gamma, Re, g_z);
        compute_RHS_3D(RHS, F, G, H, nx, ny, nz, dx, dy, dz, dt);


        poisson_3D(p, p_new, RHS, nx, ny, nz, dx, dy, dz, omega, eps, iter_max);
        inflow_condition_3D(u, v , w, nx, ny, nz, 'w', u_init, v_init, w_init);
        outflow_condition_3D(u, v, w , nx, ny, nz, 'e');
        noslip_condition_3D(u ,v , w, nx, ny, nz, 'n');
        noslip_condition_3D(u ,v , w, nx, ny, nz, 's');
        noslip_condition_3D(u ,v , w, nx, ny, nz, 'f');
        noslip_condition_3D(u ,v , w, nx, ny, nz, 'b');

        compute_uv_3D(u, v, w, F, G, H, p_new, nx, ny, nz, dx, dy, dz,dt);
//        compute_phi_3D(phi,phi_new,u,v,nx,ny,dx,dy,dt,gamma,Re) ;
//        update_3D(phi,phi_new,nx,ny) ;
//        phi_condition_3D(phi,nx,ny,dx,dy,'w',phi_init) ;
//        phi_condition_3D(phi,nx,ny,dx,dy,'e',phi_init) ;
//        phi_condition_3D(phi,nx,ny,dx,dy,'n',phi_init) ;
//        phi_condition_3D(phi,nx,ny,dx,dy,'s',phi_init) ;
        if (num_iter == 1 || num_iter%5 == 0) {
            cout << "time step : " << num_iter <<"\n" ;
            paraview_3D(num_iter ,"p" , p, nx, ny, nz, dx, dy, dz);
            paraview_3D(num_iter ,"F" , F, nx, ny, nz, dx, dy, dz);
            paraview_3D(num_iter ,"u" , u, nx, ny, nz, dx, dy, dz);
            paraview_3D(num_iter ,"v" , v, nx, ny, nz, dx, dy, dz) ;
//            paraview_3D(num_iter ,"phi" , phi, nx, ny, dx, dy);
        }

    }
}
