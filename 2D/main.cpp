#include <iostream>
#include <cmath>
using namespace std ;

void initialize(double **var, int nx, int ny , double k) ;
void visualize(double **var, int nx, int ny) ;

void update(double **var, double **var_new, int nx, int ny) ;

void paraview_vector(int num_iter, double **var_u,double **var_v,double**var_p ,double**var_phi, int nx, int ny, double dx, double dy) ;
void save_restartfile(int num_iter, const string& varName, double **var, int nx, int ny) ;
void read_restartfile(int start_num, const string& varName, double **var, int nx, int ny) ;

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
void compute_phi(double **var_phi,double **var_phi_new,double **var_u, double **var_v,  int nx, int ny,double dx,double dy,double dt,double gamma, double Re) ;


int main() {
    const int start_iter = 0 ;
    const int iter = 5000;
    const int iter_max = 1000;
    const int nx = 101; // 300
    const int ny = 11; // 10
    const double dy = 0.1;
    const double dx = 3.*dy ;
    const double dt = 0.01;
    const double u_init = 1.;
    const double v_init = 0.;
    const double p_init = 0.;
    const double phi_init = 1. ;
    const double Re = 300.;
    const double g_x = 0.;
    const double g_y = 0.;
    const double gamma = 0.5;
    const double omega = 1.7; //  SOR
    const double eps = 0.00001; // Stopping tolerance for pressure iteration

    // Construct all var 2D array
    double **u; // u^(n)
    u = (double **) malloc ((nx + 2) * sizeof(double));
    for (int i = 0; i <= nx + 1; i++) {
        u[i] = (double *) malloc ((ny + 2) * sizeof(double));
    }

    double **v; // v^(n)
    v = (double **) malloc ((nx + 2) * sizeof(double));
    for (int i = 0; i <= nx + 1; i++) {
        v[i] = (double *) malloc ((ny + 2) * sizeof(double));
    }

    double **p; // P^(n)
    p = (double **) malloc ((nx + 2) * sizeof(double));
    for (int i = 0; i <= nx + 1; i++) {
        p[i] = (double *) malloc ((ny + 2) * sizeof(double));
    }

    double **F; // F^(n)
    F = (double **) malloc ((nx + 2) * sizeof(double));
    for (int i = 0; i <= nx + 1; i++) {
        F[i] = (double *) malloc ((ny + 2) * sizeof(double));
    }

    double **G; // G^(n)
    G = (double **) malloc ((nx + 2) * sizeof(double));
    for (int i = 0; i <= nx + 1; i++) {
        G[i] = (double *) malloc ((ny + 2) * sizeof(double));
    }

    double **RHS; // RHS^(n)
    RHS = (double **) malloc ((nx + 2) * sizeof(double));
    for (int i = 0; i <= nx + 1; i++) {
        RHS[i] = (double *) malloc ((ny + 2) * sizeof(double));
    }

    double **p_new; // P^(n+1)
    p_new = (double **) malloc ((nx + 2) * sizeof(double));
    for (int i = 0; i <= nx + 1; i++) {
        p_new[i] = (double *) malloc ((ny + 2) * sizeof(double));
    }

    double **phi;
    phi = (double **) malloc ((nx + 2) * sizeof(double));
    for (int i = 0; i <= nx + 1; i++) {
        phi[i] = (double *) malloc ((ny + 2) * sizeof(double));
    }

    double **phi_new;
    phi_new = (double **) malloc ((nx + 2) * sizeof(double));
    for (int i = 0; i <= nx + 1; i++) {
        phi_new[i] = (double *) malloc ((ny + 2) * sizeof(double));
    }

    initialize(u,nx, ny,u_init) ;
    initialize(v,nx, ny,v_init);
    initialize(p,nx, ny,p_init);
    initialize(phi,nx,ny,0) ;

    // West
    inflow_condition(u, v, nx, ny, 'w', u_init, v_init);
    pressure_condition(p, nx, ny, dx, dy, 'w', 'N', p_init);
    phi_condition(phi,nx,ny,dx,dy,'w',phi_init) ;
    // East
    outflow_condition(u, v, nx, ny, 'e');
    pressure_condition(p, nx, ny, dx, dy, 'e', 'D', p_init);
    phi_condition(phi,nx,ny,dx,dy,'e',phi_init) ;
    // Top - Wall
    noslip_condition(u, v, nx, ny, 'n');
    pressure_condition(p, nx, ny, dx, dy, 'n', 'N', p_init);
    phi_condition(phi,nx,ny,dx,dy,'n',phi_init) ;
    // Bottom - Wall
    noslip_condition(u, v, nx, ny, 's');
    pressure_condition(p, nx, ny, dx, dy, 's', 'N', p_init);
    phi_condition(phi,nx,ny,dx,dy,'s',phi_init) ;
    // read file (.dat) at specific time step
    if (start_iter != 0){
        read_restartfile(start_iter,"F", F, nx, ny) ;
        read_restartfile(start_iter,"G", G, nx, ny) ;
        read_restartfile(start_iter,"RHS", RHS, nx, ny) ;
        read_restartfile(start_iter,"p", p, nx, ny) ;
        read_restartfile(start_iter,"p_new", p_new, nx, ny) ;
        read_restartfile(start_iter,"u", u, nx, ny) ;
        read_restartfile(start_iter,"v", v, nx, ny) ;
        read_restartfile(start_iter,"phi", phi, nx, ny) ;
        read_restartfile(start_iter,"phi_new", phi_new, nx, ny) ;
    }
    for (int num_iter = start_iter; num_iter <= iter; num_iter++) {

        compute_F(F, u, v, nx, ny, dx, dy, dt, gamma, Re, g_x);
        compute_G(G, u, v, nx, ny, dx, dy, dt, gamma,Re, g_y);
        compute_RHS(RHS, F, G, nx, ny, dx, dy, dt);
        // compute P at time step n+1
        poisson(p, p_new, RHS, nx, ny, dx, dy, omega, eps, iter_max);
        // compute u and v  at time step n+1
        compute_uv(u, v, F, G, p, nx, ny, dx, dy, dt);
        // apply boundary condition to u and v
        inflow_condition(u ,v , nx, ny, 'w', u_init, v_init);
        outflow_condition(u ,v , nx, ny, 'e');
        noslip_condition(u ,v , nx, ny, 'n');
        noslip_condition(u ,v , nx, ny, 's');
        // compute phi at time step n+1
        compute_phi(phi,phi_new,u,v,nx,ny,dx,dy,dt,gamma,Re) ;
        update(phi,phi_new,nx,ny) ;
        // apply boundary condition to phi
        phi_condition(phi,nx,ny,dx,dy,'w',phi_init) ;
        phi_condition(phi,nx,ny,dx,dy,'e',phi_init) ;
        phi_condition(phi,nx,ny,dx,dy,'n',phi_init) ;
        phi_condition(phi,nx,ny,dx,dy,'s',phi_init) ;

        if(num_iter%1000 == 0){
            // save file (.dat) every 1000 time step
            save_restartfile(num_iter,"F", F, nx, ny) ;
            save_restartfile(num_iter,"G", G, nx, ny) ;
            save_restartfile(num_iter,"RHS", RHS, nx, ny) ;
            save_restartfile(num_iter,"p", p, nx, ny) ;
            save_restartfile(num_iter,"p_new", p_new, nx, ny) ;
            save_restartfile(num_iter,"u", u, nx, ny) ;
            save_restartfile(num_iter,"v", v, nx, ny) ;
            save_restartfile(num_iter,"phi", phi, nx, ny) ;
            save_restartfile(num_iter,"phi_new", phi_new, nx, ny) ;

        }
        if (num_iter == 0 || num_iter%5 == 0) {
            cout << "time step : " << num_iter <<"\n" ;
//            visualize(u,nx,ny);
            paraview_vector(num_iter, u, v, p, phi, nx, ny, dx, dy) ;
        }
    }
}
