#include <iostream>
#include <cmath>
using namespace std ;

void initialize_3D(double ***var, int nx, int ny ,int nz, double c) ;
void visualize_3D(double ***var, int nx, int ny, int nz) ;
void update_3D(double ***var, double ***var_new, int nx, int ny, int nz) ;
void paraview_3D(int num_iter, const string& varName, double ***var, int nx, int ny, int nz, double dx, double dy,double dz) ;
void paraview_vector_3D(int num_iter, double ***var_u,double ***var_v,double ***var_w,double***var_p ,double***var_phi, int nx, int ny,int nz, double dx, double dy, double dz) ;
void save_restartfile_3D(int num_iter, const string& varName, double ***var, int nx, int ny, int nz) ;
void read_restartfile_3D(int start_num, const string& varName, double ***var, int nx, int ny, int nz) ;

void noslip_condition_3D(double ***var_u, double ***var_v,double ***var_w, int nx, int ny, int nz, char side);
void outflow_condition_3D(double ***var_u, double ***var_v,double ***var_w, int nx, int ny, int nz , char side) ;
void inflow_condition_3D(double ***var_u, double ***var_v,double ***var_w, int nx, int ny, int nz, char side, double u, double v, double w);
void pressure_condition_3D(double ***var_P, int nx, int ny, int nz, char side) ;
void phi_condition_3D(double ***var_phi,int nx,int ny, int nz, double D_jett, double phi) ;

void compute_F_3D(double ***var_F, double ***var_u, double ***var_v, double ***var_w, int nx, int ny, int nz, double dx, double dy, double dz, double dt, double gamma, double Re, double g_x);
void compute_G_3D(double ***var_G, double ***var_u, double ***var_v, double ***var_w, int nx, int ny, int nz, double dx, double dy, double dz, double dt, double gamma, double Re, double g_y);
void compute_H_3D(double ***var_H, double ***var_u, double ***var_v, double ***var_w, int nx, int ny, int nz, double dx, double dy, double dz, double dt, double gamma, double Re, double g_z);
void compute_RHS_3D(double ***var_RHS,double ***var_F, double ***var_G, double ***var_H, int nx, int ny, int nz, double dx, double dy, double dz, double dt);

void poisson_3D(double ***var_p, double ***var_p_new, double ***RHS, int nx, int ny, int nz, double dx, double dy, double dz, double omega, double eps, int iter_max);


void compute_uv_3D(double ***var_u, double ***var_v, double ***var_w, double ***var_F, double ***var_G, double ***var_H, double ***var_p_new, int nx, int ny, int nz, double dx, double dy, double dz, double dt) ;
void compute_phi_3D(double ***var_phi,double ***var_phinew,double ***var_u, double ***var_v, double ***var_w, int nx, int ny, int nz, double dx, double dy, double dz, double dt,double gamma, double Re) ;


int main() {
    const int start_iter = 0 ;
    const int iter = 1;
    const int iter_max = 1000;
    const int D  = 1;
    const int nx = 20 * D;
    const int ny = 15 * D;
    const int nz = 10 * D;
    const double dy = 1;
    const double dx = dy ;
    const double dz = dy ;// 3
    const double dt = 0.1;
    const double u_init = 0.5;
    const double v_init = 0.;
    const double w_init = 0.;
    const double p_init = 0.;
    const double phi_init = 1. ;
    const double v_jett = 1. ;
    const double D_jett = 3*D ;
    const double Re = 1000.; // Reynolds number
    const double g_x = 0.;
    const double g_y = 0.;
    const double g_z = 0.;
    const double gamma = 0.5; // Upwind differencing factor
    const double omega = 1.7; // Relaxation parameter for SOR iteration
    const double eps = 0.00001; // Stopping tolerance for pressure iteration
    int i , j  ;

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
    initialize_3D(u,nx,ny,nz,u_init) ;
    initialize_3D(v,nx,ny,nz,0) ;
    initialize_3D(w,nx,ny,nz,0) ;
    initialize_3D(p,nx,ny,nz,0) ;
    initialize_3D(phi,nx,ny,nz,0) ;
    // Left
    inflow_condition_3D(u, v, w, nx, ny, nz, 'w', u_init, v_init, w_init);
    pressure_condition_3D(p, nx, ny, nz, 'w');
    // Right
    outflow_condition_3D(u, v, w, nx, ny, nz, 'e');
    pressure_condition_3D(p, nx, ny, nz, 'e');
    // Top
    outflow_condition_3D(u, v, w, nx, ny, nz, 'n');
    pressure_condition_3D(p, nx, ny, nz, 'n');
    // South
    noslip_condition_3D(u,v,w, nx, ny, nz, 's');
    pressure_condition_3D(p, nx, ny, nz, 's');
    phi_condition_3D(phi, nx, ny, nz, D_jett, phi_init) ;
    // Front
    outflow_condition_3D(u, v,w, nx, ny, nz, 'f');
    pressure_condition_3D(p, nx, ny, nz, 'f');
    // Behind
    outflow_condition_3D(u, v,w, nx, ny, nz, 'b');
    pressure_condition_3D(p, nx, ny, nz, 'b');

//    if (start_iter != 0){
//        read_restartfile_3D(start_iter,"F", F, nx, ny, nz) ;
//        read_restartfile_3D(start_iter,"G", G, nx, ny, nz) ;
//        read_restartfile_3D(start_iter,"H", H, nx, ny, nz) ;
//        read_restartfile_3D(start_iter,"RHS", RHS, nx, ny, nz) ;
//        read_restartfile_3D(start_iter,"p", p, nx, ny, nz) ;
//        read_restartfile_3D(start_iter,"p_new", p_new, nx, ny, nz) ;
//        read_restartfile_3D(start_iter,"u", u, nx, ny, nz) ;
//        read_restartfile_3D(start_iter,"v", v, nx, ny, nz) ;
//        read_restartfile_3D(start_iter,"w", w, nx, ny, nz) ;
//        read_restartfile_3D(start_iter,"phi", phi, nx, ny, nz) ;
//        read_restartfile_3D(start_iter,"phi_new", phi_new, nx, ny, nz) ;
//    }
    for (int num_iter = start_iter + 1 ; num_iter <= iter; num_iter++) {
        cout << "time step : " << num_iter <<"\n" ;
        compute_F_3D(F, u, v, w, nx, ny, nz, dx, dy, dz, dt, gamma, Re, g_x);
        compute_G_3D(G, u, v, w, nx, ny, nz, dx, dy, dz, dt, gamma, Re, g_y);
        compute_H_3D(H, u, v, w, nx, ny, nz, dx, dy, dz, dt, gamma, Re, g_z);
        compute_RHS_3D(RHS, F, G, H, nx, ny, nz, dx, dy, dz, dt);
        cout  << " ----------------------u-----------------------------\n" ;
        for( j = ny+1 ; j >= 0 ; j--){
            for(i = 0 ; i <= nx + 1 ; i++){
                cout << u[i][j][nz/2] << " " ;
            }
            cout << "\n" ;
        }
        cout  << " -----------------------P----------------------------\n" ;
        for( j = ny+1 ; j >= 0 ; j--){
            for(i = 0 ; i <= nx + 1 ; i++){
                cout << p[i][j][nz/2] << " " ;
            }
            cout << "\n" ;
        }
        cout  << " ----------------------P_new-----------------------------\n" ;
        for( j = ny+1 ; j >= 0 ; j--){
            for(i = 0 ; i <= nx + 1 ; i++){
                cout << p_new[i][j][nz/2] << " " ;
            }
            cout << "\n" ;
        }
        poisson_3D(p, p_new, RHS, nx, ny, nz, dx, dy, dz, omega, eps, iter_max);
        inflow_condition_3D(u, v , w, nx, ny, nz, 'w', u_init, v_init, w_init);
        outflow_condition_3D(u, v, w , nx, ny, nz, 'e');
        outflow_condition_3D(u ,v , w, nx, ny, nz, 'n');
        noslip_condition_3D(u ,v , w, nx, ny, nz, 's');
        outflow_condition_3D(u ,v , w, nx, ny, nz, 'f');
        outflow_condition_3D(u ,v , w, nx, ny, nz, 'b');
        compute_uv_3D(u, v, w, F, G, H, p_new, nx, ny, nz, dx, dy, dz,dt);
        compute_phi_3D(phi,phi_new,u,v,w,nx,ny,nz,dx,dy,nz,dt,gamma,Re) ;
        update_3D(phi, phi_new, nx, ny, nz) ;
        phi_condition_3D(phi, nx, ny, nz, D_jett, phi_init) ;

        paraview_vector_3D(num_iter, u, v, w, p, phi, nx, ny, nz, dx, dy, dz) ;
        if (num_iter == 1 || num_iter%5 == 0) {
            save_restartfile_3D(num_iter,"F", F, nx, ny, nz) ;
            save_restartfile_3D(num_iter,"G", G, nx, ny, nz) ;
            save_restartfile_3D(num_iter,"H", H, nx, ny, nz) ;
            save_restartfile_3D(num_iter,"RHS", RHS, nx, ny, nz) ;
            save_restartfile_3D(num_iter,"p", p, nx, ny, nz) ;
            save_restartfile_3D(num_iter,"p_new", p_new, nx, ny, nz) ;
            save_restartfile_3D(num_iter,"u", u, nx, ny, nz) ;
            save_restartfile_3D(num_iter,"v", v, nx, ny, nz) ;
            save_restartfile_3D(num_iter,"w", w, nx, ny, nz) ;
            save_restartfile_3D(num_iter,"phi", phi, nx, ny, nz) ;
            save_restartfile_3D(num_iter,"phi_new", phi_new, nx, ny, nz) ;
        }



    }
}
