#include <iostream>
#include <cmath>
using namespace std ;

void initialize(double **var, int nx, int ny) ;
void visualize(double **var, int nx, int ny) ;
void update(double **var, double **var_new, int nx, int ny) ;
void paraview(int i, string varName, double **var, int nx, int ny, int dx, int dy);

//library_boundary.cpp
void noslip_condition(double **var_u, double **var_v, int nx, int ny, char wall);
void outflow_condition(double **var_u, double **var_v, int nx, int ny, char wall);
void inflow_condition(double **var_u, double **var_v, int nx, int ny, char wall, double u, double v);
void pressure_condition(double **var_P, int nx, int ny, double dx, double dy, char wall, char type, double P) ;
void phi_condition(double **var_phi, int nx, int ny) ;

//library_comp.cpp
void compute_F(double **var_F, double **var_u, double **var_v, int nx, int ny, double dx, double dy, double dt, double gamma, double Re, double g_x);
void compute_G(double **var_G, double **var_u, double **var_v, int nx, int ny, double dx, double dy, double dt, double gamma, double Re, double g_y);
void compute_RHS(double **var_RHS,double **var_F, double **var_G, int nx, int ny, double dx, double dy, double dt) ;

//library_poisson.cpp
void poisson(double **var_p, double **var_p_new, double **RHS, int nx, int ny, double dx, double dy, double omega, double eps, int iter_max);

//library_simulation.cpp
void compute_uv(double **var_u, double **var_v, double **var_F, double **var_G, double **var_p_new, int nx, int ny, double dx, double dy, double dt);



int main() {
    const int iter(5000);
    const int iter_max(1000);
    const int nx(400);
    const int ny(20);
    const double dx = 1.;
    const double dy = 1.;
    const double dt = 0.01;
    const double u_ini = 1.;
    const double v_ini = 0.;
    const double p_ini = 0.;

    const double Re = 1.; // Reynolds number
    const double g_x = 0.;
    const double g_y = 0.;
    const double gamma = 0.5; // Upwind differencing factor
    const double omega = 1.7; // Relaxation parameter for SOR iteration
    const double esp = 0.00001; // Stopping tolerance for pressure iteration

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

    double **phi_new; // P^(n+1)
    phi_new = (double **) malloc ((nx + 2) * sizeof(double));
    for (int i = 0; i <= nx + 1; i++) {
        phi_new[i] = (double *) malloc ((ny + 2) * sizeof(double));
    }

    initialize(u,nx, ny) ;
    initialize(v,nx, ny);
    initialize(p,nx, ny);

    // West
    inflow_condition(u, v, nx, ny, 'w', u_ini, v_ini);
    pressure_condition(p, nx, ny, dx, dy, 'w', 'n', 0);
    // East
    outflow_condition(u, v, nx, ny, 'e');
    pressure_condition(p, nx, ny, dx, dy, 'e', 'd', 0);
    // Wall (North and South)
    noslip_condition(u, v, nx, ny, 'n');
    pressure_condition(p, nx, ny, dx, dy, 'n', 'n', 0);
    noslip_condition(u, v, nx, ny, 's');
    pressure_condition(p, nx, ny, dx, dy, 's', 'n', 0);
    // cout << " - Init BC - " << endl;
    phi_condition(phi, nx,ny) ;
    for (int i = 1; i <= iter; i++) {

        compute_F(F, u, v, nx, ny, dx, dy, dt, gamma, Re, g_x);
        compute_G(G, u, v, nx, ny, dx, dy, dt, gamma,Re, g_y);
        compute_RHS(RHS, F, G, nx, ny, dx, dy, dt);

        poisson(p, p_new, RHS, nx, ny, dx, dy, omega, esp, iter_max);
        inflow_condition(u ,v , nx, ny, 'w', 1., 0.);
        outflow_condition(u ,v , nx, ny, 'e');
        noslip_condition(u ,v , nx, ny, 'n');
        noslip_condition(u ,v , nx, ny, 's');

        compute_uv(u, v, F, G, p_new, nx, ny, dx, dy, dt);
        visualize(u, nx, ny) ;
        if (i == 1 || i%5 == 0) {
            cout << "Iter "<< i << endl;
            paraview(i ,"p" , p, nx, ny, dx, dy);
            paraview(i ,"u" , u, nx, ny, dx, dy);
            paraview(i ,"v" , u, nx, ny, dx, dy) ;
            paraview(i ,"phi" , phi, nx, ny, dx, dy);
        }

    }
}
