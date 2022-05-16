#include <cmath>
#include <iostream>
using namespace std ;
int main (){
    const double pi = atan(1)*4 ;
    const int D = 4 ;
    const int grid_size_x = 20*D  ;
    const int grid_size_y = 10*D  ;
    const int x_c = 5*D ;
    const int y_c = 5*D ;
    const int radius = 5 ;
    int d_i ;
    double x, y, d, angle, p_w, p_h, a ;
    int phi[grid_size_x+1][grid_size_y+1] ;
    for (int i = 0; i <= grid_size_x; i++ ){
        for(int j = 0 ; j <= grid_size_y; j++){
            phi[i][j] = 0 ;
        }
    }
    for (int col=0; col <= grid_size_x;col++) {
        for (int row=0; row <= grid_size_y; row++) {
            x = x_c - col;
            y = y_c - row;
            d = sqrt(pow(x,2) + pow(y,2)); // Distance from center
            d_i = floor(d); // Integer part of distance
            if (d_i < radius){
                phi[col][row] = 1 ;
            } else if (d_i > radius){
                phi[col][row] = 0;
            } else{
                // Pixel is touching circumference line.
                // Lets calculate how much area is inside circle.
                // Get angle from center.
                if (x == 0) {
                    if (y > 0) {
                        angle = pi;
                    } else {
                        angle = -pi;
                    }
                } else {
                    angle = atan(y / x);
                }
                p_w = (d - d_i) * sin(angle);
                p_h = (d - d_i) * cos(angle);
                a = p_w * p_h; // (rough) Pixel area under circumference, which falls inside circle.
                if (a >= 0.50) { phi[col][row] = 1 ;}
                else {phi[col][row] = 0 ;}
            }
        }
    }
    for (int j = 0; j <= grid_size_y; j++ ){
        for(int i = 0 ; i <= grid_size_x; i++){
            cout << phi[i][j] << " " ;
        }
        cout << "\n" ;
    }
}



