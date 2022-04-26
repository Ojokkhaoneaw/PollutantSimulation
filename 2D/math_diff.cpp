#include <cmath>
using namespace std;

double du_dx(double **var_u, int i, int j, double dx) {
  return (var_u[i][j]-var_u[i-1][j])/dx;
}

double du2_dx(double **var_u, int i, int j, double dx, double gamma) {
  double Lterm = pow((var_u[i][j]+var_u[i+1][j])/2, 2.) - pow((var_u[i-1][j]+var_u[i][j])/2, 2.) ;
  double Rterm = 0.25*gamma*( abs(var_u[i][j]+var_u[i+1][j])*(var_u[i][j]-var_u[i+1][j]) - abs(var_u[i-1][j]+var_u[i][j])*(var_u[i-1][j]-var_u[i][j]) );
  return (Rterm+Lterm)/dx;
}

double duv_dy(double **var_u, double **var_v, int i, int j, double dy, double gamma) {
  double Lterm = 0.25*( (var_v[i][j]+var_v[i+1][j])*(var_u[i][j]+var_u[i][j+1]) - (var_v[i][j-1]+var_v[i+1][j-1])*(var_u[i][j-1]+var_u[i][j]) );
  double Rterm = 0.25*gamma*( abs(var_v[i][j]+var_v[i+1][j])*(var_u[i][j]-var_u[i][j+1]) - abs(var_v[i][j-1]+var_v[i+1][j-1])*(var_u[i][j-1]-var_u[i][j]) );
  return (Rterm+Lterm)/dy;
}

double d2u_dx2(double **var_u, int i, int j, double dx) {
  return (var_u[i+1][j]-(2*var_u[i][j])+var_u[i-1][j])/(dx*dx);
}

double d2u_dy2(double **var_u, int i, int j, double dy) {
  return (var_u[i][j+1]-2*var_u[i][j]+var_u[i][j-1])/(dy*dy);
}

double dp_dx(double **var_p, int i, int j, double dx) {
  return (var_p[i+1][j]-var_p[i][j])/dx;
}

double dv_dy(double **var_v, int i, int j, double dy) {
  return (var_v[i][j]-var_v[i][j-1])/dy;
}

double duv_dx(double **var_u, double **var_v, int i, int j, double dx, double gamma) {
  double Lterm = 0.25*( (var_u[i][j]+var_u[i][j+1])*(var_v[i][j]+var_v[i+1][j]) - (var_u[i-1][j]+var_u[i-1][j+1])*(var_v[i-1][j]+var_v[i][j]) );
  double Rterm = 0.25*gamma*( abs(var_u[i][j]+var_u[i][j+1])*(var_v[i][j]-var_v[i+1][j]) - abs(var_u[i-1][j]+var_u[i-1][j+1])*(var_v[i-1][j]-var_v[i][j]) );
  return (Lterm+Rterm)/dx;
}

double dv2_dy(double **var_v, int i, int j, double dy, double gamma) {
  double Lterm =  pow((var_v[i][j]+var_v[i][j+1])/2, 2.) - pow((var_v[i][j-1]+var_v[i][j])/2, 2.) ;
  double Rterm = 0.25*gamma*( abs(var_v[i][j]+var_v[i][j+1])*(var_v[i][j]-var_v[i][j+1]) - abs(var_v[i][j-1]+var_v[i][j])*(var_v[i][j-1]-var_v[i][j]) );
  return (Lterm+Rterm)/dy;
}

double d2v_dx2(double **var_v, int i, int j, double dx) {
  return (var_v[i+1][j]-2*var_v[i][j]+var_v[i-1][j])/(dx*dx);
}

double d2v_dy2(double **var_v, int i, int j, double dy) {
  return (var_v[i][j+1]-2*var_v[i][j]+var_v[i][j-1])/(dy*dy);
}

double dp_dy(double **var_p, int i, int j, double dy) {
  return (var_p[i][j+1]-var_p[i][j])/dy;
}

double dphi_dx(double **var_phi, int i, int j, double dx){
  return (var_phi[i+1][j]-var_phi[i][j])/dx;  
}


double dphi_dy(double **var_phi, int i, int j, double dy){
  return (var_phi[i][j+1]-var_phi[i][j])/dy;  
}

double d2phi_dx2(double **var_phi, int i, int j, double dx){
  return (var_phi[i+1][j]-2*var_phi[i][j]+var_phi[i-1][j])/(dx*dx);
}

double d2phi_dy2(double **var_phi, int i, int j, double dy){
  return (var_phi[i][j+1]-2*var_phi[i][j]+var_phi[i][j-1])/(dy*dy);
}

double dvphi_dy(double **var_v, double **var_phi, int i, int j, double dy, double gamma) {
  double Lterm = 0.5*( var_v[i][j]*(var_phi[i][j]+var_phi[i][j+1]) -  var_v[i][j-1]*(var_phi[i][j-1]+var_phi[i][j]) );
  double Rterm = 0.5*gamma*( abs(var_v[i][j])*(var_phi[i][j]-var_phi[i][j+1]) -  abs(var_v[i][j-1])*(var_phi[i][j-1]-var_phi[i][j]));
  return (Rterm+Lterm)/dy;
}

double duphi_dx(double **var_u, double **var_phi, int i, int j, double dx, double gamma) {
  double Lterm = 0.5*( var_u[i][j]*(var_phi[i][j]+var_phi[i+1][j]) -  var_u[i-1][j]*(var_phi[i-1][j]+var_phi[i][j]) );
  double Rterm = 0.5*gamma*( abs(var_u[i][j])*(var_phi[i][j]-var_phi[i+1][j]) -  abs(var_u[i-1][j])*(var_phi[i-1][j]-var_phi[i][j]));
  return (Rterm+Lterm)/dx;
}
