# include <iostream>
# include <fstream>
# include <cmath>
using namespace std ;

// Initialize variable
void initialize(double **var, int nx, int ny) {
    for (int i = 0 ; i <= nx-1 ; i++){
        for (int j = 0 ; j <= ny-1 ; j++) {
            var[i][j] = 0.0 ;
        }
    }
}

// Visualize variable
void visualize(double **var, int nx, int ny) {
    for (int i = 0 ; i <= nx-1 ; i++){
        for (int j = 0 ; j <= ny-1 ; j++) {
            cout << var[i][j] << " " ;
        }
        cout << "\n" ;
    }
    cout << "\n" ;
}



//Update variable
void update(double **var, double **var_new, int nx, int ny) {
    for (int i = 0 ; i <= nx-1 ; i++){
        for (int j = 0 ; j <= ny-1 ; j++) {
            var[i][j] = var_new[i][j] ;
        }
    }
}
// incomplete
void paraview(double **var,int nx, int ny) {
    ofstream myfile ;
    myfile.open("var.vtk") ;
//  Paraview header
    myfile << "# vtk DataFile Version 2.0\n" ;
    myfile << "FlowField\n" ;
    myfile << "ACSII\n" ;
//  Grid
    myfile << "DATASET STRUCTURED_GRID\n" ;
    myfile << "DIMENSIONS" << nx << " " << 1 << " " << ny << "\n" ;
    myfile << "POINTS" << nx*1*ny << "float\n" ;
    for(int j = 0 ; j <= ny-1 ; j++ ){
        for(int i = 0 ; i <= nx - 1 ; i++ ) {
            myfile << var[i][j] << " " ;
        }
        myfile << "\n" ;
    }
    myfile << "\n" ;
    myfile << "POINT_DATA" ;
    myfile << "POINTS" << nx*ny << "\n" ;
    myfile << "\n" ;
    myfile << "SCALARS PHI float 1\n" ;
    myfile << "LOOKUP_TABLE default\n" ;
    for(int j = 0 ; j <= ny-1 ; j++) {
        for(int i = 0 ; i <= nx-1 ; i++) {
            myfile << var[i][j] << "\n" ;
        }
    }
    myfile.close() ;
}



// Simulate 2D heat equation
void simulation(double **var, double **var_new, int nx, int ny, double c, double dx, double dy, double dt){
    double RHS ;
    // Simulation
    for (int i = 1 ; i <= nx-2 ; i++) {
        for(int j = 1 ; j <= ny-2; j++) {
            RHS = c * ((var[i +1][j] - 2.*var[i][j] + var[i-1][j]) / (dx* dy))+
                  c * ((var[i][j+1] - 2.*var[i][j] + var[i][j-1]) / (dx* dy)) ;
            var_new[i][j] = var[i][j] + dt * RHS;
        }
    }
}

// Boundary Condition
//1. Inflow - > Dirichlet Inflow Boundary Condition (u0 = 1 , phi(y<=0.5) = 1, phi(y>0.5) = 0)
//2. Wall - > Dirichlet & No-slip Boundary Condition (Bottom Wall) (u = v = 0)
//3. Outflow - > Neumann Outflow Boundary Condition (du_dx = dv_dx = 0)
void noslip_condition(double **var_u, double **var_v, int nx, int ny, char side) {
    switch (side) {
        case 'w': // Left Side of Grid
            for (int j = 1; j <= ny + 1; j++) {
                var_u[0][j] = 0;
                var_v[0][j] = -var_v[1][j];
            }
            break;
        case 'n': // Top Side of Grid
            for (int i = 1; i <= nx + 1; i++) {
                var_v[i][ny] = 0;
                var_u[i][ny+1] = -var_u[i][ny] ;
            }
            break;
        case 'e': // Right Side of Grid
            for (int j = 1; j <= ny + 1; j++) {
                var_v[nx+1][j] = -var_v[nx][j];
                var_u[nx][j] = 0;
            }
            break;
        case 's': // Bottom Side of Grid
            for (int i = 1; i <= nx + 1; i++) {
                var_v[i][0] = 0;
                var_u[i][0] = -var_u[i][1];
            }
            break;
    }
}

void freeslip_condition(double **var_u, double **var_v, int nx, int ny, char side) {
    switch (side) {
        case 'w':
            for (int j = 1; j <= ny + 1; j++) {
                var_u[0][j] = 0;
                var_v[0][j] = var_v[1][j];
            }
            break;
        case 'n':
            for (int i = 1; i <= nx + 1; i++) {
                var_u[i][ny+1] = var_u[i][ny];
                var_v[i][ny] = 0;
            }
            break;
        case 'e':
            for (int j = 1; j <= ny + 1; j++) {
                var_u[nx][j] = 0;
                var_v[nx+1][j] = var_v[nx][j];
            }
            break;
        case 's':
            for (int i = 1; i <= nx + 1; i++) {
                var_u[i][0] = var_u[i][1];
                var_v[i][0] = 0;
            }
            break;
    }
}

void outflow_condition(double **var_u, double **var_v, int nx, int ny, char side) {
    switch (side) {
        case 'w':
            for (int j = 1; j <= ny + 1; j++) {
                var_u[0][j] = var_u[1][j];
                var_v[0][j] = var_v[1][j];
            }
            break;
        case 'n':
            for (int i = 1; i <= nx + 1; i++) {
                var_u[i][ny+1] = var_u[i][ny];
                var_v[i][ny+1] = var_v[i][ny-1];
            }
            break;
        case 'e':
            for (int j = 1; j <= ny + 1; j++) {
                var_u[nx+1][j] = var_u[nx-1][j];
                var_v[nx+1][j] = var_v[nx][j];
            }
            break;
        case 's':
            for (int i = 1; i <= nx + 1; i++) {
                var_u[i][0] = var_u[i][1];
                var_v[i][0] = var_v[i][1];
            }
            break;
    }
}

void inflow_condition(double **var_u, double **var_v, int nx, int ny, char side, double u, double v) {
    switch (side) {
        case 'w':
            for (int j = 1; j <= ny + 1; j++) {
                var_u[0][j] = u;
                var_v[0][j] = 2*v - var_v[1][j];
            }
            break;
        case 'n':
            for (int i = 1; i <= nx + 1; i++) {
                var_u[i][ny+1] = 2*u - var_u[i][ny];
                var_v[i][ny] = v;
            }
            break;
        case 'e':
            for (int j = 1; j <= ny + 1; j++) {
                var_u[nx][j] = u;
                var_v[nx+1][j] = 2*v - var_v[nx][j];
            }
            break;
        case 's':
            for (int i = 1; i <= nx + 1; i++) {
                var_u[i][0] = 2*u - var_u[i][1];
                var_v[i][1] = v;
            }
            break;
    }
}

// void boundary_periodic(double **var_u, double **var_v, int nx, int ny, char wall) {

// }

void pressure_condition(double **var_P, int nx, int ny, double dx, double dy, char side, char type, double P) {
    switch (type) {
        case 'D': // Dirichlet Boundary Condition : First-Type
            switch (side) {
                case 'w':
                    for (int j = 1; j <= ny + 1; j++) {
                        var_P[0][j] = 2*P - var_P[1][j];
                    }
                    break;
                case 'n':
                    for (int i = 1; i <= nx + 1; i++) {
                        var_P[i][ny+1] = 2*P - var_P[i][ny];
                    }
                    break;
                case 'e':
                    for (int j = 1; j <= ny + 1; j++) {
                        var_P[nx+1][j] = 2*P - var_P[nx][j];
                    }
                    break;
                case 's':
                    for (int i = 1; i <= nx + 1; i++) {
                        var_P[i][0] = 2*P - var_P[i][1];
                    }
                    break;
            }
            break;
        case 'N': // Neumann Boundary Condition : Second-Type
            switch (side) {
                case 'w':
                    for (int j = 1; j <= ny + 1; j++) {
                        var_P[0][j] = var_P[1][j] - dx*P;
                    }
                    break;
                case 'n':
                    for (int i = 1; i <= nx + 1; i++) {
                        var_P[i][0] = var_P[i][1] + dx*P;
                    }
                    break;
                case 'e':
                    for (int j = 1; j <= ny + 1; j++) {
                        var_P[nx+1][j] = var_P[nx][j] + dx*P;
                    }
                    break;
                case 's':
                    for (int i = 1; i <= nx + 1; i++) {
                        var_P[i][ny+1] = var_P[i][ny] - dx*P;
                    }
                    break;
            }
            break;
    }
}

void phi_boundary(double **var_phi, int nx, int ny) {
    for(int j = 0; j <= ny+1 ; j++) {
        if (j <= 0.5) {
            var_phi[0][j] = 1 ;
        }else{
            var_phi[0][j] = 0 ;
        }
    }
}

