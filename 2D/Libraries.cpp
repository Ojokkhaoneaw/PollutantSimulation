# include <iostream>
# include <fstream>
# include <cmath>
using namespace std ;

// Initialize variable
void initialize(double **var, int nx, int ny , double c) {
    for (int i = 0 ; i <= nx+1 ; i++){
        for (int j = 0 ; j <= ny+1 ; j++) {
            var[i][j] = c ;
        }
    }
}
// Initialize variable
void initialize3D(double ***var, int nx, int ny ,int nz, double c) {
    for (int i = 0 ; i <= nx+1 ; i++){
        for (int j = 0 ; j <= ny+1 ; j++) {
            for (int k = 0 ; k <= nz+1 ; k++) {
                var[i][j][k] = c ;
            }
        }
    }
}

// Visualize variable
void visualize(double **var, int nx, int ny) {
    for (int j = ny + 1; j >= 0; j--){
        for (int i = 0 ; i <= nx+1 ; i++) {
            cout << var[i][j] << " " ;
        }
        cout << "\n" ;
    }
    cout << "\n" ;
}

void visualize3D(double ***var, int nx, int ny, int nz){
    for (int i = 0; i < nx+1; ++i) {
        for (int j = 0; j < ny+1; ++j) {
            for (int k = 0; k < nz+1; ++k) {
                cout << var[i][j][k] << " " ;
            }
            cout << "\n" ;
        }
        cout << "\n" ;
    }
    cout << "\n" ;
}



//Update variable
void update(double **var, double **var_new, int nx, int ny) {
    for (int i = 0 ; i <= nx+1 ; i++){
        for (int j = 0 ; j <= ny+1 ; j++) {
            var[i][j] = var_new[i][j] ;
        }
    }
}
// incomplete

void paraview(int num_iter, const string& varName, double **var, int nx, int ny, double dx, double dy) {
    string fileName = "var_" + varName + "_" + to_string(num_iter) + ".vtk";
    ofstream myfile;
    myfile.open(fileName);

    //Paraview Header
    myfile << "# vtk DataFile Version 2.0" << endl;
    myfile << "FlowField" << endl;
    myfile << "ASCII" << endl;

    //Grid
    myfile << "DATASET STRUCTURED_GRID" << endl;
    myfile << "DIMENSIONS"<<" "<< nx+2 << " " << ny+2 << " " << 1 << " " << endl;
    myfile << "POINTS"<<" "<< (nx+2)*(ny+2) <<" " << "double" << endl;
    for (int j = 0; j <= ny + 1; j++) {
        for (int i = 0; i <= nx + 1; i++) {
            myfile << i*dx << " " << j*dy <<" "<< "0" << endl;
        }
    }
    myfile << endl;

    //Data
    myfile << "POINT_DATA"<< " " << (nx+2)*(ny+2) << endl;
    myfile << endl;
    myfile << "SCALARS"<< " " << varName <<" " <<"double" << endl;
    myfile << "LOOKUP_TABLE default" << endl;
    for (int j = 0; j <= ny+1; j++) {
        for (int i = 0; i <= nx+1; i++) {
            myfile << var[i][j] << endl;
        }
    }
}

void paraview3D(int num_iter, const string& varName, double ***var, int nx, int ny, int nz, double dx, double dy,double dz) {
    string fileName = "var3D_" + varName + "_" + to_string(num_iter) + ".vtk";
    ofstream myfile;
    myfile.open(fileName);

    //Paraview Header
    myfile << "# vtk DataFile Version 2.0" << endl;
    myfile << "FlowField" << endl;
    myfile << "ASCII" << endl;

    //Grid
    myfile << "DATASET STRUCTURED_GRID" << endl;
    myfile << "DIMENSIONS"<<" "<< nx+2 << " " << ny+2 << " " << nz+2 << " " << endl;
    myfile << "POINTS"<<" "<< (nx+2)*(ny+2)*(nz+2) <<" " << "double" << endl;
    for (int i = 0; i <= nx+1; i++) {
        for (int j = 0; j <= ny+1; j++) {
            for(int k = 0; k <= nz+1; k++){
                myfile << i*dx << " " << j*dy <<" "<< k*dz << endl;
            }
        }
    }
    myfile << endl;

    //Data
    myfile << "POINT_DATA"<< " " << (nx+2)*(ny+2)*(nz+2) << endl;
    myfile << endl;
    myfile << "SCALARS"<< " " << varName <<" " <<"double" << endl;
    myfile << "LOOKUP_TABLE default" << endl;
    for (int i = 0; i <= nx+1; i++) {
        for (int j = 0; j <= ny+1; j++) {
            for(int k = 0; k <= nz+1; k++){
                myfile << var[i][j][k] << endl;
            }
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
            for (int j = 0; j <= ny + 1; j++) {
                var_v[0][j] = -var_v[1][j];
            }
            break;
        case 'e': // Right Side of Grid
            for (int j = 0; j <= ny + 1; j++) {
                var_v[nx+1][j] = -var_v[nx][j];
            }
            break;
        case 'n': // Top Side of Grid
            for (int i = 0; i <= nx + 1; i++) {
                var_u[i][ny+1] = -var_u[i][ny] ;
            }
            break;
        case 's': // Bottom Side of Grid
            for (int i = 0; i <= nx + 1; i++) {
                var_u[i][0] = -var_u[i][1] ;
            }
            break;
    }
}


void outflow_condition(double **var_u, double **var_v, int nx, int ny, char side) {
    switch (side) {
        case 'w':
            for (int j = 0; j <= ny + 1; j++) {
                var_u[0][j] = var_u[1][j];
                var_v[0][j] = var_v[1][j];
            }
            break;
        case 'e':
            for (int j = 0; j <= ny + 1; j++) {
                var_u[nx][j] = var_u[nx-1][j];
                var_v[nx+1][j] = var_v[nx][j];
            }
            break;
        case 'n':
            for (int i = 0; i <= nx + 1; i++) {
                var_u[i][ny+1] = var_u[i][ny];
                var_v[i][ny] = var_v[i][ny-1];
            }
            break;
        case 's':
            for (int i = 0; i <= nx + 1; i++) {
                var_u[i][0] = var_u[i][1];
                var_v[i][0] = var_v[i][1];
            }
            break;
    }
}

void inflow_condition(double **var_u, double **var_v, int nx, int ny, char side, double u, double v) {
    switch (side) {
        case 'w':
            for (int j = 0; j <= ny + 1; j++) {
                var_u[0][j] = u;
                var_v[0][j] = 2*v - var_v[1][j];
            }
            break;
        case 'e':
            for (int j = 0; j <= ny + 1; j++) {
                var_u[nx][j] = u;
                var_v[nx+1][j] = 2*v - var_v[nx][j];
            }
            break;
        case 'n':
            for (int i = 0; i <= nx + 1; i++) {
                var_u[i][ny+1] = 2*u - var_u[i][ny];
                var_v[i][ny] = v;
            }
            break;
        case 's':
            for (int i = 0; i <= nx + 1; i++) {
                var_u[i][0] = 2*u - var_u[i][1];
                var_v[i][1] = v;
            }
            break;
    }
}


void pressure_condition(double **var_P, int nx, int ny, double dx, double dy, char side, char type, double P) {
    switch (type) {
        case 'D': // Dirichlet Boundary Condition : First-Type
            switch (side) {
                case 'w':
                    for (int j = 0; j <= ny + 1; j++) {
                        var_P[0][j] = 2*P - var_P[1][j];
                    }
                    break;
                case 'n':
                    for (int i = 0; i <= nx + 1; i++) {
                        var_P[i][ny+1] = 2*P - var_P[i][ny];
                    }
                    break;
                case 'e':
                    for (int j = 0; j <= ny + 1; j++) {
                        var_P[nx+1][j] = 2*P - var_P[nx][j];
                    }
                    break;
                case 's':
                    for (int i = 0; i <= nx + 1; i++) {
                        var_P[i][0] = 2*P - var_P[i][1];
                    }
                    break;
            }
            break;
        case 'N': // Neumann Boundary Condition : Second-Type
            switch (side) {
                case 'w':
                    for (int j = 0; j <= ny + 1; j++) {
                        var_P[0][j] = var_P[1][j];
                    }
                    break;
                case 'n':
                    for (int i = 0; i <= nx + 1; i++) {
                        var_P[i][0] = var_P[i][1];
                    }
                    break;
                case 'e':
                    for (int j = 0; j <= ny + 1; j++) {
                        var_P[nx+1][j] = var_P[nx][j];
                    }
                    break;
                case 's':
                    for (int i = 0; i <= nx + 1; i++) {
                        var_P[i][ny+1] = var_P[i][ny];
                    }
                    break;
            }
            break;
    }
}

void phi_condition(double **var_phi,int nx,int ny,double dx,double dy, char side, double phi) {
    switch (side) {
        case 'w':
            double phi_c ;
            for (int j = 0; j <= ny + 1; j++) {
                if (j <= (ny+1)/2 ) {phi_c = phi ;}
                else{ phi_c = 0;}
                var_phi[0][j] = 2*phi_c - var_phi[1][j];
            }
            break;
        case 'n':
            for (int i = 0; i <= nx + 1; i++) {
                var_phi[i][0] = var_phi[i][1];
            }
            break;
        case 'e':
            for (int j = 0; j <= ny + 1; j++) {
                var_phi[nx+1][j] = var_phi[nx][j];
            }
            break;
        case 's':
            for (int i = 0; i <= nx + 1; i++) {
                var_phi[i][ny+1] = var_phi[i][ny];
            }
            break;
    }
}

