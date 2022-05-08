#include <cmath>
using namespace std;

//-------------------x-axis---------------

double du_dx(double ***var_u, int i, int j, int k,double dx) {
  return (var_u[i][j][k]-var_u[i-1][j][k])/dx;
}

double du2_dx(double ***var_u, int i, int j, int k,double dx, double gamma) {
  double Lterm = pow((var_u[i][j][k]+var_u[i+1][j][k])/2, 2.) - pow((var_u[i-1][j][k]+var_u[i][j][k])/2, 2.) ;
  double Rterm = 0.25*gamma*( abs(var_u[i][j][k]+var_u[i+1][j][k])*(var_u[i][j][k]-var_u[i+1][j][k]) - abs(var_u[i-1][j][k]+var_u[i][j][k])*(var_u[i-1][j][k]-var_u[i][j][k]) );
  return (Rterm+Lterm)/dx;
}

double duv_dy(double ***var_u, double ***var_v, int i, int j, int k, double dy, double gamma) {
  double Lterm = 0.25*( (var_v[i][j][k]+var_v[i+1][j][k])*(var_u[i][j][k]+var_u[i][j+1][k]) - (var_v[i][j-1][k]+var_v[i+1][j-1][k])*(var_u[i][j-1][k]+var_u[i][j][k]) );
  double Rterm = 0.25*gamma*( abs(var_v[i][j][k]+var_v[i+1][j][k])*(var_u[i][j][k]-var_u[i][j+1][k]) - abs(var_v[i][j-1][k]+var_v[i+1][j-1][k])*(var_u[i][j-1][k]-var_u[i][j][k]) );
  return (Rterm+Lterm)/dy;
}

double duw_dz(double ***var_u, double ***var_w, int i, int j, int k, double dz, double gamma) {
  double Lterm = 0.25*( (var_w[i][j][k]+var_w[i+1][j][k])*(var_u[i][j][k]+var_u[i][j][k+1]) - (var_w[i][j][k-1]+var_w[i+1][j][k-1])*(var_u[i][j][k-1]+var_u[i][j][k]) );
  double Rterm = 0.25*gamma*( abs(var_w[i][j][k]+var_w[i+1][j][k])*(var_u[i][j][k]-var_u[i][j][k+1]) - abs(var_w[i][j][k-1]+var_w[i+1][j][k-1])*(var_u[i][j][k-1]-var_u[i][j][k]) );
  return (Rterm+Lterm)/dz;
}

double d2u_dx2(double ***var_u, int i, int j, int k, double dx) {
  return (var_u[i+1][j][k]-(2*var_u[i][j][k])+var_u[i-1][j][k])/(dx*dx);
}

double d2u_dy2(double ***var_u, int i, int j, int k, double dy) {
  return (var_u[i][j+1][k]-2*var_u[i][j][k]+var_u[i][j-1][k])/(dy*dy);
}


double d2u_dz2(double ***var_u, int i, int j, int k, double dz) {
  return (var_u[i][j][k+1]-2*var_u[i][j][k]+var_u[i][j][k-1])/(dz*dz);

}
double dp_dx(double ***var_p, int i, int j, int k, double dx) {
  return (var_p[i+1][j][k]-var_p[i][j][k])/dx;
}

//------------------- y -axis--------------
double dv_dy(double ***var_v, int i, int j, int k, double dy) {
  return (var_v[i][j][k]-var_v[i][j-1][k])/dy;
}

double duv_dx(double ***var_u, double ***var_v, int i, int j, int k, double dx, double gamma) {
  double Lterm = 0.25*( (var_u[i][j][k]+var_u[i][j+1][k])*(var_v[i][j][k]+var_v[i+1][j][k]) - (var_u[i-1][j][k]+var_u[i-1][j+1][k])*(var_v[i-1][j][k]+var_v[i][j][k]) );
  double Rterm = 0.25*gamma*( abs(var_u[i][j][k]+var_u[i][j+1][k])*(var_v[i][j][k]-var_v[i+1][j][k]) - abs(var_u[i-1][j][k]+var_u[i-1][j+1][k])*(var_v[i-1][j][k]-var_v[i][j][k]) );
  return (Lterm+Rterm)/dx;
}

double dvw_dz(double ***var_v, double ***var_w, int i, int j, int k, double dz, double gamma) {
  double Lterm = 0.25*( (var_v[i][j][k]+var_v[i][j][k+1])*(var_w[i][j][k]+var_w[i][j+1][k]) - (var_v[i][j-1][k]+var_v[i][j-1][k+1])*(var_w[i][j-1][k]+var_w[i][j][k]) );
  double Rterm = 0.25*gamma*( abs(var_v[i][j][k]+var_v[i][j][k+1])*(var_w[i][j][k]-var_w[i][j+1][k]) - abs(var_v[i][j-1][k]+var_v[i][j-1][k+1])*(var_w[i][j-1][k]-var_w[i][j][k]) );
  return (Lterm+Rterm)/dz;
}

double dv2_dy(double ***var_v, int i, int j, int k, double dy, double gamma) {
  double Lterm =  pow((var_v[i][j][k]+var_v[i][j+1][k])/2, 2.) - pow((var_v[i][j-1][k]+var_v[i][j][k])/2, 2.) ;
  double Rterm = 0.25*gamma*( abs(var_v[i][j][k]+var_v[i][j+1][k])*(var_v[i][j][k]-var_v[i][j+1][k]) - abs(var_v[i][j-1][k]+var_v[i][j][k])*(var_v[i][j-1][k]-var_v[i][j][k]) );
  return (Lterm+Rterm)/dy;
}

double d2v_dx2(double ***var_v, int i, int j, int k, double dx) {
  return (var_v[i+1][j][k]-2*var_v[i][j][k]+var_v[i-1][j][k])/(dx*dx);
}

double d2v_dy2(double ***var_v, int i, int j, int k, double dy) {
  return (var_v[i][j+1][k]-2*var_v[i][j][k]+var_v[i][j-1][k])/(dy*dy);
}

double d2v_dz2(double ***var_v, int i, int j, int k, double dz) {
  return (var_v[i][j][k+1]-2*var_v[i][j][k]+var_v[i][j][k-1])/(dz*dz);
}

double dp_dy(double ***var_p, int i, int j, int k, double dy) {
  return (var_p[i][j+1][k]-var_p[i][j][k])/dy;
}

//---------------- z-axis ------------------------
double dw_dz(double ***var_w, int i, int j, int k, double dz) {
  return (var_w[i][j][k]-var_w[i][j][k-1])/dz;
}

double dw2_dz(double ***var_w, int i, int j, int k, double dz, double gamma) {
  double Lterm = pow((var_w[i][j][k]+var_w[i][j][k+1])/2, 2.) - pow((var_w[i][j][k-1]+var_w[i][j][k])/2, 2.) ;
  double Rterm = 0.25*gamma*( abs(var_w[i][j][k]+var_w[i][j][k+1])*(var_w[i][j][k]-var_w[i][j][k+1]) - abs(var_w[i][j][k-1]+var_w[i][j][k])*(var_w[i][j][k-1]-var_w[i][j][k]) );
  return (Rterm+Lterm)/dz;
}

double duw_dx(double ***var_u, double ***var_w, int i, int j, int k, double dx, double gamma) {
  double Lterm = 0.25*( (var_u[i][j][k]+var_u[i][j][k+1])*(var_w[i][j][k]+var_w[i+1][j][k]) - (var_u[i-1][j][k]+var_u[i-1][j][k+1])*(var_w[i-1][j][k]+var_w[i][j][k]) );
  double Rterm = 0.25*gamma*( abs(var_u[i][j][k]+var_u[i][j][k+1])*(var_w[i][j][k]-var_w[i+1][j][k]) - abs(var_u[i-1][j][k]+var_u[i-1][j][k+1])*(var_w[i-1][j][k]-var_w[i][j][k]) );
  return (Rterm+Lterm)/dx;
}

double dvw_dy(double ***var_v, double ***var_w, int i, int j, int k, double dx, double gamma) {
  double Lterm = 0.25*( (var_v[i][j][k]+var_v[i][j][k+1])*(var_w[i][j][k]+var_w[i][j+1][k]) - (var_v[i][j-1][k]+var_v[i][j-1][k+1])*(var_w[i][j-1][k]+var_w[i][j][k]) );
  double Rterm = 0.25*gamma*( abs(var_v[i][j][k]+var_v[i][j][k+1])*(var_w[i][j][k]-var_w[i][j+1][k]) - abs(var_v[i][j-1][k]+var_v[i][j-1][k+1])*(var_w[i][j-1][k]-var_w[i][j][k]) );
  return (Rterm+Lterm)/dy;
}

double d2w_dx2(double ***var_w, int i, int j,int k, double dx) {
  return (var_w[i+1][j][k]-(2*var_w[i][j][k])+var_w[i-1][j][k])/(dx*dx);
}

double d2w_dy2(double ***var_w, int i, int j, int k, double dy) {
  return (var_w[i][j+1][k]-(2*var_w[i][j][k])+var_w[i][j-1][k])/(dy*dy);
}

double d2w_dz2(double ***var_w, int i, int j, int k, double dz) {
  return (var_w[i][j][k+1]-2*var_w[i][j][k]+var_w[i][j][k-1])/(dz*dz);
}

double dp_dz(double ***var_p, int i, int j, int k, double dz) {
  return (var_p[i][j][k+1]-var_p[i][j][k])/dz;
}


double dphi_dx(double ***var_phi, int i, int j, int k, double dx){
  return (var_phi[i+1][j][k]-var_phi[i][j][k])/dx;  
}


double dphi_dy(double ***var_phi, int i, int j, int k, double dy){
  return (var_phi[i][j+1][k]-var_phi[i][j][k])/dy;  
}


double dphi_dz(double ***var_phi, int i, int j, int k, double dz){
  return (var_phi[i][j][k+1]-var_phi[i][j][k])/dz;  
}

double d2phi_dx2(double ***var_phi, int i, int j, int k, double dx){
  return (var_phi[i+1][j][k]-2*var_phi[i][j][k]+var_phi[i-1][j][k])/(dx*dx);
}

double d2phi_dy2(double ***var_phi, int i, int j, int k, double dy){
  return (var_phi[i][j+1][k]-2*var_phi[i][j][k]+var_phi[i][j-1][k])/(dy*dy);
}


double d2phi_dz2(double ***var_phi, int i, int j, int k, double dz){
  return (var_phi[i][j][k+1]-2*var_phi[i][j][k]+var_phi[i][j][k-1])/(dz*dz);
}


double dvphi_dy(double ***var_v, double ***var_phi, int i, int j, int k, double dy, double gamma) {
  double Lterm = 0.5*( var_v[i][j][k]*(var_phi[i][j][k]+var_phi[i][j+1][k]) -  var_v[i][j-1][k]*(var_phi[i][j-1][k]+var_phi[i][j][k]) );
  double Rterm = 0.5*gamma*( abs(var_v[i][j][k])*(var_phi[i][j][k]-var_phi[i][j+1][k]) -  abs(var_v[i][j-1][k])*(var_phi[i][j-1][k]-var_phi[i][j][k]));
  return (Rterm+Lterm)/dy;
}

double duphi_dx(double ***var_u, double ***var_phi, int i, int j, int k, double dx, double gamma) {
  double Lterm = 0.5*( var_u[i][j][k]*(var_phi[i][j][k]+var_phi[i+1][j][k]) -  var_u[i-1][j][k]*(var_phi[i-1][j][k]+var_phi[i][j][k]) );
  double Rterm = 0.5*gamma*( abs(var_u[i][j][k])*(var_phi[i][j][k]-var_phi[i+1][j][k]) -  abs(var_u[i-1][j][k])*(var_phi[i-1][j][k]-var_phi[i][j][k]));
  return (Rterm+Lterm)/dx;
}


double dwphi_dz(double ***var_w, double ***var_phi, int i, int j, int k, double dx, double gamma) {
  double Lterm = 0.5*( var_w[i][j][k]*(var_phi[i][j][k]+var_phi[i][j][k+1]) -  var_w[i][j][k-1]*(var_phi[i][j][k-1]+var_phi[i][j][k]) );
  double Rterm = 0.5*gamma*( abs(var_w[i][j][k])*(var_phi[i][j][k]-var_phi[i][j][k+1]) -  abs(var_w[i][j][k-1])*(var_phi[i][j][k-1]-var_phi[i][j][k]));
  return (Rterm+Lterm)/dz;
}
