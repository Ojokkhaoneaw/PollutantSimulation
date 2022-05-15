# include <iostream>
# include <fstream>
# include <cmath>
using namespace std ;

// Initialize variable
void initialize_3D(double ***var, int nx, int ny ,int nz, double c) {
    for (int i = 0 ; i <= nx+1 ; i++){
        for (int j = 0 ; j <= ny+1 ; j++) {
            for (int k = 0 ; k <= nz+1 ; k++) {
                var[i][j][k] = c ;
            }
        }
    }
}

//Update variable
void update_3D(double ***var, double ***var_new, int nx, int ny, int nz) {
    for (int i = 0 ; i <= nx+1 ; i++){
        for (int j = 0 ; j <= ny+1 ; j++){
            for(int k = 0 ; k <= nz+1 ; k++) {
                var[i][j][k] = var_new[i][j][k] ;
            }
        }
    }
}

void paraview_3D(int num_iter, const string& varName, double ***var, int nx, int ny, int nz, double dx, double dy,double dz) {
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
    for (int k = 0; k <= nz+1; k++) {
        for (int j = 0; j <= ny+1; j++) {
            for(int i = 0; i<= nx+1; i++){
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
    for (int k = 0; k <= nz+1; k++) {
        for (int j = 0; j <= ny+1; j++) {
            for(int i = 0; i<= nx+1; i++){
                myfile << var[i][j][k] << endl;
            }
        }
    }
    myfile.close() ;
}

void save_restartfile_3D(int num_iter, const string& varName, double ***var, int nx, int ny, int nz){
    string fileName = "var3D_" + varName + "_" + to_string(num_iter) + ".dat" ;
    ofstream myfileO; // output file stream
    myfileO.open(fileName);
    for(int k = 0 ; k <= nz+1 ; k++){
        for(int j = 0; j <= ny+1; j++){
            for(int i = 0; i <= nx+1; i++){
                myfileO << var[i][j][k] << endl;
            }
        }
    }
    myfileO.close();
}

void read_restartfile_3D(int start_num, const string& varName, double ***var, int nx, int ny, int nz){
    string fileName = "var3D_" + varName + "_" + to_string(start_num) + ".dat" ;
    ifstream myfileI; // input file stream
    myfileI.open(fileName);
    for(int k = 0 ; k <= nz+1 ; k++){
        for(int j = 0; j <= ny+1; j++){
            for(int i = 0; i <= nx+1; i++) {
                myfileI >> var[i][j][k];
            }
        }
    }
    myfileI.close();
}



// Boundary Condition
//1. Inflow - > Dirichlet Inflow Boundary Condition (u0 = 1 , phi(y<=0.5) = 1, phi(y>0.5) = 0)
//2. Wall - > Dirichlet & No-slip Boundary Condition (Bottom Wall) (u = v = 0)
//3. Outflow - > Neumann Outflow Boundary Condition (du_dx_3D = dv_dx = 0)

void noslip_condition_3D(double ***var_u, double ***var_v,double ***var_w, int nx, int ny, int nz, char side) {
    switch (side) {
        case 'w': // yz plane of grid cell
            for (int j = 1; j <= ny ; j++) {
                for(int k = 1 ; k <= nz ; k++ ) {
                    var_u[0][j][k] = 0 ;
                    var_v[0][j][k] = -var_v[1][j][k];
                    var_w[0][j][k] = -var_w[1][j][k];
                }
            }
            break;
        case 'e': //
            for (int j = 1; j <= ny ; j++) {
                for(int k = 1 ; k <= nz ; k++ ) {
                    var_u[nx + 1][j][k] = 0 ;
                    var_v[nx + 1][j][k] = -var_v[nx][j][k];
                    var_w[nx + 1][j][k] = -var_w[nx][j][k];
                }
            }
            break;
        case 'n': // xz plane of grid cell
            for (int i = 1; i <= nx; i++) {
                for (int k = 1; k <= nz; k++) {
                    var_u[i][ny + 1][k] = -var_u[i][ny][k];
                    var_v[i][ny + 1][k] = 0;
                    var_w[i][ny + 1][k] = -var_w[i][ny][k]  ;
                }
            }
            break;
        case 's': //
            for (int i = 1; i <= nx; i++) {
                for(int k = 1 ; k <= nz ; k++ ) {
                    var_u[i][0][k] = -var_u[i][1][k];
                    var_v[i][0][k] = 0;
                    var_w[i][0][k] = -var_w[i][1][k] ;
                }
            }
            break;
        case 'f': // xy plane of grid cell
            for (int i = 1; i <= nx; i++) {
                for(int j = 1 ; j <= ny ; j++ ) {
                    var_u[i][j][0] = -var_u[i][j][1];
                    var_v[i][j][0] = -var_w[i][j][1];
                    var_w[i][j][0] = 0;
                }
            }
            break;
        case 'b': //
            for (int i = 1; i <= nx; i++) {
                for(int j = 1 ; j <= ny ; j++ ) {
                    var_u[i][j][nz+1] = -var_u[i][j][nz];
                    var_w[i][j][nz+1] = -var_w[i][j][nz];
                    var_w[i][j][nz+1] = 0;
                }
            }
            break;
    }
}

void outflow_condition_3D(double ***var_u, double ***var_v,double ***var_w, int nx, int ny, int nz , char side) {
    switch (side) {
        case 'w': // yz plane of grid cell
            for (int j = 1; j <= ny; j++) {
                for(int k = 1 ; k <= nz ; k++ ) {
                    var_u[0][j][k] = var_u[1][j][k];
                    var_v[0][j][k] = var_v[1][j][k];
                    var_w[0][j][k] = var_w[1][j][k];
                }
            }
            break;
        case 'e':
            for (int j = 1; j <= ny; j++) {
                for(int k = 1 ; k <= nz ; k++ ) {
                    var_u[nx][j][k]   = var_u[nx - 1][j][k];
                    var_v[nx+1][j][k] = var_v[nx][j][k];
                    var_w[nx+1][j][k] = var_w[nx][j][k];
                }
            }
            break;
        case 'n': // xz plane of grid cell
            for (int i = 1; i <= nx; i++) {
                for(int k = 1 ; k <= ny ; k++ ) {
                    var_u[i][ny+1][k]  = var_u[i][ny][k];
                    var_v[i][ny][k]    = var_v[i][ny-1][k];
                    var_w[i][ny+1][k]  = var_w[i][ny][k];
                }
            }
            break;
        case 's':
            for (int i = 1; i <= nx; i++) {
                for(int k = 1 ; k <= nz ; k++ ) {
                    var_u[i][0][k] = var_u[i][1][k];
                    var_v[i][0][k] = var_v[i][1][k];
                    var_w[i][0][k] = var_w[i][1][k];
                }
            }
            break;
        case 'f': // xy plane of grid cell
            for (int i = 1; i <= nx; i++) {
                for(int j = 1 ; j <= ny ; j++ ) {
                    var_u[i][j][0] = var_u[i][j][1];
                    var_v[i][j][0] = var_v[i][j][1];
                    var_w[i][j][0] = var_w[i][j][1];
                }
            }
            break;
        case 'b': // Back Side of Grid
            for (int i = 1; i <= nx; i++) {
                for(int j = 1 ; j <= ny ; j++ ){
                    var_u[i][j][nz+1] = var_u[i][j][nz];
                    var_v[i][j][nz+1] = var_v[i][j][nz];
                    var_w[i][j][nz]   = var_w[i][j][nz-1];
                }
            }
            break;
    }
}

void inflow_condition_3D(double ***var_u, double ***var_v,double ***var_w, int nx, int ny, int nz, char side, double u, double v, double w) {
    switch (side) {
        case 'w':
            for (int j = 1; j <= ny; j++) {
                for(int k = 1 ; k <= nz ; k++ ) {
                    var_u[0][j][k] = u;
                    var_v[0][j][k] = 2 * v - var_v[1][j][k];
                    var_w[0][j][k] = 2 * w - var_w[1][j][k];
                }
            }
            break;
        case 'e':
            for (int j = 1; j <= ny; j++) {
                for(int k = 1 ; k <= nz ; k++ ) {
                    var_u[nx][j][k] = u;
                    var_v[nx + 1][j][k] = 2 * v - var_v[nx][j][k];
                    var_w[nx + 1][j][k] = 2 * w - var_w[nx][j][k];
                }
            }
            break;
        case 'n':
            for (int i = 1; i <= nx; i++) {
                for(int k = 1 ; k <= nz ; k++ ) {
                    var_u[i][ny + 1][k] = 2 * u - var_u[i][ny][k];
                    var_v[i][ny][k] = v;
                    var_w[i][ny + 1][k] = 2 * w - var_w[i][ny][k];
                }
            }
            break;
        case 's':
            for (int i = 1; i <= nx; i++) {
                for (int k = 1; k <= nz; k++) {
                    var_u[i][0][k] = 2 * u - var_u[i][1][k];
                    var_v[i][0][k] = v;
                    var_w[i][0][k] = 2 * w - var_w[i][1][k];
                }
            }
            break;
        case 'f':
            for (int i = 1; i <= nx; i++) {
                for(int j = 1 ; j <= ny ; j++ ) {
                    var_u[i][j][nz+1] = 2 * u - var_u[i][j][nz];
                    var_v[i][j][nz+1] = 2 * v - var_v[i][j][nz];
                    var_w[i][j][nz] = w;
                }
            }
            break;
        case 'b':
            for (int i = 1; i <= nx; i++) {
                for (int j = 1; j <= ny; j++) {
                    var_u[i][j][0] = 2 * u - var_u[i][j][1];
                    var_v[i][j][0] = 2 * v - var_v[i][j][1];
                    var_w[i][j][0] = w;
                }
            }
            break;
    }
}


void pressure_condition_3D(double ***var_P, int nx, int ny, int nz,  double dx, double dy, double dz, char side, char type, double P) {
    switch (type) {
        case 'D': // Dirichlet Boundary Condition : First-Type
            switch (side) {
                case 'w':
                    for (int j = 1; j <= ny; j++) {
                        for (int k = 1; k <= nz; k++) {
                            var_P[0][j][k] = 2 * P - var_P[1][j][k];
                        }
                    }
                    break;
                case 'n':
                    for (int i = 1; i <= nx; i++) {
                        for (int k = 1; k <= nz; k++) {
                            var_P[i][ny + 1][k] = 2 * P - var_P[i][ny][k];
                        }
                    }
                    break;
                case 'e':
                    for (int j = 1; j <= ny; j++) {
                        for (int k = 1; k <= nz; k++) {
                            var_P[nx + 1][j][k] = 2 * P - var_P[nx][j][k];
                        }
                    }
                    break;
                case 's':
                    for (int i = 1; i <= nx; i++) {
                        for (int k = 1; k <= nz; k++) {
                            var_P[i][0][k] = 2 * P - var_P[i][1][k];
                        }
                    }
                    break;
                case 'f':
                    for (int i = 1; i <= nx; i++) {
                        for (int j = 1; j <= ny; j++) {
                            var_P[i][j][0] = 2 * P - var_P[i][j][1];
                        }
                    }
                    break;
                case 'b':
                    for (int i = 1; i <= nx; i++) {
                        for (int j = 1; j <= ny; j++) {
                            var_P[i][j][nz + 1] = 2 * P - var_P[i][j][nz];
                        }
                    }
                    break;
            }
            break;
        case 'N': // Neumann Boundary Condition : Second-Type
            switch (side) {
                case 'w':
                    for (int j = 1; j <= ny; j++) {
                        for (int k = 1; k <= nz; k++) {
                            var_P[0][j][k] = var_P[1][j][k];
                        }
                    }
                    break;
                case 'e':
                    for (int j = 1; j <= ny; j++) {
                        for (int k = 1; k <= nz; k++) {
                            var_P[nx + 1][j][k] = var_P[nx][j][k];
                        }
                    }
                    break;
                case 'n':
                    for (int i = 1; i <= nx; i++) {
                        for (int k = 1; k <= nz; k++) {
                            var_P[i][ny+1][k] = var_P[i][ny][k];
                        }
                    }
                    break;
                case 's':
                    for (int i = 1; i <= nx; i++) {
                        for (int k = 1; k <= nz; k++) {
                            var_P[i][0][k] = var_P[i][1][k];
                        }
                    }
                    break;
                case 'f':
                    for (int i = 1; i <= nx; i++) {
                        for (int j = 1; j <= ny; j++) {
                            var_P[i][j][nz+1] = var_P[i][j][nz];
                        }
                    }
                    break;
                case 'b':
                    for (int i = 1; i <= nx; i++) {
                        for (int j = 0; j <= ny; j++) {
                            var_P[i][j][0] = var_P[i][j][1];
                        }
                    }
                    break;
            }
        }
}


//void phi_condition_3D(double ***var_phi,int nx,int ny,double dx,double dy, char side, double phi) {
//    switch (side) {
//        case 'w':
//            double phi_c ;
//            for (int j = 0; j <= ny + 1; j++) {
//                if (j <= (ny+1)/2 ) {phi_c = phi ;}
//                else{ phi_c = 0;}
//                var_phi[0][j] = 2*phi_c - var_phi[1][j];
//            }
//            break;
//        case 'n':
//            for (int i = 0; i <= nx + 1; i++) {
//                var_phi[i][0] = var_phi[i][1];
//            }
//            break;
//        case 'e':
//            for (int j = 0; j <= ny + 1; j++) {
//                var_phi[nx+1][j] = var_phi[nx][j];
//            }
//            break;
//        case 's':
//            for (int i = 0; i <= nx + 1; i++) {
//                var_phi[i][ny+1] = var_phi[i][ny];
//            }
//            break;
//    }
//}

