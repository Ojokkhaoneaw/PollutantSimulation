# include <iostream>
# include <fstream>
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
        cout << "\n" ; ;
    }
}

//Update variable
void update(double **var, double **var_new, int nx, int ny) {
    for (int i = 0 ; i <= nx+1 ; i++){
        for (int j = 0 ; j <= ny+1 ; j++) {
            var[i][j] = var_new[i][j] ;
        }
    }
}

void paraview(double **var,int nx, int ny, double dx, double dy) {
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
            myfile << dx*i << " " << dy*j << " 0\n" ;
        }
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