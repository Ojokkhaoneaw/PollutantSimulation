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


//Update variable
void update(double **var, double **var_new, int nx, int ny) {
    for (int i = 0 ; i <= nx+1 ; i++){
        for (int j = 0 ; j <= ny+1 ; j++) {
            var[i][j] = var_new[i][j] ;
        }
    }
}


void paraview_vector(int num_iter, double **var_u,double **var_v,double**var_p ,double**var_phi, int nx, int ny, double dx, double dy) {
    string fileName = "result_" + to_string(num_iter) + ".vtk";
    ofstream myfile;
    myfile.open(fileName);

    //Paraview Header
    myfile << "# vtk DataFile Version 3.0" << endl;
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
    myfile << "SCALARS"<< " " << "pressure" <<" " <<"double" << endl;
    myfile << "LOOKUP_TABLE" << " " << "pressure" << endl;
    for (int j = 0; j <= ny+1; j++) {
        for (int i = 0; i <= nx+1; i++) {
            myfile << var_p[i][j] << endl;
        }
    }
    myfile << "SCALARS"<< " " << "phi" <<" " <<"double" << endl;
    myfile << "LOOKUP_TABLE" << " " << "phi" << endl;
    for (int j = 0; j <= ny+1; j++) {
        for (int i = 0; i <= nx+1; i++) {
            myfile << var_phi[i][j] << endl;
        }
    }
    myfile << "VECTORS"<< " " << "velocity" <<" " <<"double" << endl;
    for (int j = 0; j <= ny+1; j++) {
        for (int i = 0; i <= nx+1; i++) {
            myfile << var_u[i][j] << " " <<var_v[i][j] << " " <<"0"<< endl;
        }
    }
    myfile.close();
}
void save_restartfile(int num_iter, const string& varName, double **var, int nx, int ny){
    string fileName = "var_" + varName + "_" + to_string(num_iter) + ".dat" ;
    ofstream myfileO; // output file stream
    myfileO.open(fileName);
    for (int j = 0; j <= ny+1; j++){
        for (int i = 0 ; i <= nx+1 ; i++) {
            myfileO << var[i][j] << " " ;
        }
        myfileO << "\n" ;
    }
    myfileO << "\n" ;
    myfileO.close();
}

void read_restartfile(int start_num, const string& varName, double **var, int nx, int ny){
    string fileName = "var_" + varName + "_" + to_string(start_num) + ".dat" ;
    ifstream myfileI; // input file stream
    myfileI.open(fileName);
    for (int j = ny + 1; j >= 0; j--){
        for (int i = 0 ; i <= nx+1 ; i++) {
            myfileI >> var[i][j] ;
        }
    }
    myfileI.close();
}

// Boundary Condition
//1. Inflow - > Dirichlet Inflow Boundary Condition (u0 = 1 , phi(y<=0.5) = 1, phi(y>0.5) = 0)
//2. Wall - > Dirichlet & No-slip Boundary Condition (Bottom Wall) (u = v = 0)
//3. Outflow - > Neumann Outflow Boundary Condition (du_dx = dv_dx = 0)
void noslip_condition(double **var_u, double **var_v, int nx, int ny, char side) {
    switch (side) {
        case 'w': // Left Side of Grid
            for (int j = 0; j <= ny + 1; j++) {
                var_u[0][j] = 0;
                var_v[0][j] = -var_v[1][j];
            }
            break;
        case 'e': // Right Side of Grid
            for (int j = 0; j <= ny + 1; j++) {
                var_u[nx][j] = 0;
                var_v[nx+1][j] = -var_v[nx][j];
            }
            break;
        case 'n': // Top Side of Grid
            for (int i = 0; i <= nx + 1; i++) {
                var_u[i][ny+1] = -var_u[i][ny] ;
                var_v[i][ny] = 0;
            }
            break;
        case 's': // Bottom Side of Grid
            for (int i = 0; i <= nx + 1; i++) {
                var_u[i][0] = -var_u[i][1] ;
                var_v[i][0] = 0;
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
                var_v[i][0] = v;
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
                        var_P[i][ny+1] = var_P[i][ny];
                    }
                    break;
                case 'e':
                    for (int j = 0; j <= ny + 1; j++) {
                        var_P[nx+1][j] = var_P[nx][j];
                    }
                    break;
                case 's':
                    for (int i = 0; i <= nx + 1; i++) {
                        var_P[i][0] = var_P[i][1];
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

