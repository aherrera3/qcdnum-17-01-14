/*
Script for the DGLAP evolution of the proton PDFs using qcdnum software.

Initial PDFs are given at the value q = 2.56 GeV^2.

Proton PDFs parametrizations taken from: Bonvini, M., & Giuli, F. (2019). A new simple PDF parametrization: improved description
																										 of the HERA data. The European Physical Journal Plus, 134(10), 531.
*/

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include "QCDNUM/QCDNUM.h"
using namespace std;

// x Parton Distribution functions:
//----------------------------------------------------------------------
double xupv(double x)
{
  double A_uv = 10.19304899633023;
  double B_uv = 0.76, C_uv = 4.6, E_uv = 2.6, F_uv = 0.35, G_uv = 0.049;
  double pd = A_uv * pow(x, B_uv) * pow((1 - x), C_uv) * (1 + E_uv * pow(x, 2) + F_uv * log(x) + G_uv * pow(log(x), 2));
  return pd;
}
//----------------------------------------------------------------------
double xdnv(double x)
{
  double A_dv = 5.564255181349489, B_dv = 0.99, C_dv = 4.7;
  double pd = A_dv * pow(x, B_dv) * pow((1 - x), C_dv);
  return pd;
}
//----------------------------------------------------------------------
double xglu(double x)
{
  double A_g = 0.872978687751462, B_g = -0.52, C_g = 4.5, F_g = 0.217, G_g = 0.0112;
  double pd = A_g * pow(x, B_g) * pow((1 - x), C_g) * (1 + F_g * log(x) + G_g * pow(log(x), 2));
  return pd;
}
//----------------------------------------------------------------------
double xdbar(double x)
{
  double A_dbar = 0.14, B_dbar = -0.33, C_dbar = 24, D_dbar = 38, F_dbar = 0.071;
  double pd = A_dbar * pow(x, B_dbar) * pow((1 - x), C_dbar) * (1 + D_dbar * x + F_dbar * log(x));
  return pd;
}
//----------------------------------------------------------------------
double xubar(double x)
{
  double A_ubar = 0.14, B_ubar = -0.33, C_ubar = 11, D_ubar = 18, F_ubar = 0.071;
  double pd = A_ubar * pow(x, B_ubar) * pow((1 - x), C_ubar) * (1 + D_ubar * x + F_ubar * log(x));
  return pd;
}
//----------------------------------------------------------------------
double xsbar(double x)
{
  double f_s = 0.4;
  double gamma_s = f_s / (1 - f_s);
  double pd = gamma_s * xdbar(x);
  return pd;
}

//----------------------------------------------------------------------
/*
** Function that calls all xPDF(x) at initial scale q0.
*/
double func(int *ipdf, double *x)
{
  int i = *ipdf;
  double xb = *x;
  double f = 0;
  if (i == 0)
    f = xglu(xb);   // gluon
  if (i == 1)
    f = xdnv(xb);   // valence down
  if (i == 2)
    f = xupv(xb);   // valence up
  if (i == 3)
    f = 0;          // valence strange
  if (i == 4)
    f = xdbar(xb);  // Sea quarks:
  if (i == 5)
    f = xubar(xb);
  if (i == 6)
    f = xsbar(xb);
  if (i == 7)
    f = 0;        // charm
  if (i == 8)
    f = 0;        // anti-charm
  if (i == 9)
    f = 0;        // bottom
  if (i == 10)
    f = 0;        // anti-bottom
  if (i == 11)
    f = 0;        // top
  if (i == 12)
    f = 0;        // anti-top
  return f;
}

int main()
{
  // to delete the old .csv files
  system("rm /opt/qcdnum-17-01-14/output/*.csv");

  // unpolarized pdfset, NNLO, VFNS
  int pdfset_type = 1, iord = 3, nfix = 0;

  // Generating x and q arrays logarithmic
  const int n_points = 200;
  double x[n_points];

  double e_mi = -3, e_ma = 0;
  double delta_e = (e_ma - e_mi)/(n_points-1);
  int i = 0;
  for (double e=e_mi; e<e_ma; e+=delta_e){
    x[i] = pow(10, e);
    //cout << x[i] << endl;
    i++;
  }
  
  // q^2 array
  static const double q2i = 0.5, q2f = 10, delta_q2 = 0.1;
  static const int n_q2 =  int((q2f-q2i)/delta_q2 +0.5);
  cout << "n_q2= " <<  n_q2 << endl;
  double q2[n_q2];

  for (int i=0; i<n_q2; i++){
    q2[i] = pow(2.56, 2.0*(q2i+i*delta_q2));
    //cout << q2[i] << endl;
  }

  //double q2[] = {2.56, 5, 10, 100, 140, 170, 1e3, 5e3, 7e3,
  //               1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10,// 1e11, 1e12};


  // size of x and y arrays, and min(x)
  int nx = size(x);
  double q20 = 0.5, q2max = q2[n_q2 - 1];    // initial and final boundaries of q^2 qcdnum grid

  // Quarks flavour composition: is an input for evolfg:
  double pdf_flavour[] =
      // tb  bb  cb  sb  ub  db   g   d   u   s   c   b   t
        {0., 0., 0., 0., 0., -1., 0., 1., 0., 0., 0., 0., 0., // 1=dval    // because dv = d - dbar
        0., 0., 0., 0., -1., 0., 0., 0., 1., 0., 0., 0., 0., // 2=uval
        0., 0., 0., -1., 0., 0., 0., 0., 0., 1., 0., 0., 0., // 3=sval
        0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0.,  // 4=dbar
        0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0.,  // 5=ubar
        0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0.,  // 6=sbar
        0., 0., -1., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., // 7=cval
        0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,  // 8=cbar
        0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,  // 9=zero
        0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,  // 10=zero
        0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,  // 11=zero
        0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}; // 12=zero

  // variables for the qcinit function
  int lun = -6;
  string outfile = " "; // if lun=6, qcdnum messages appear in the standard output

  // variables for the gxmake function (x-grid parameters) 
  //double xmin[] = x;
  int iwt[] = {1,  1,  1,  1, 1};            // x grid weigths
  int n = size(iwt);                      // = size(xmin)             
  int nxout;                      
  int nxin = 200, iosp = 2;           // nxin: number of x grid points. 
                                      // iosp: order for spline interpolation (2=linear, 3=quadratic). 

  // variables for the gqmake function (q^2-grid parameters)
  double qarr[] = {q20, q2max};
  //cout << q20 << " " << q2max << endl;
  double wgt[] = {1., 1.};   
  int nqarr = size(qarr);            // = size(wgt)      
  int nqin = 120, nqout;             //  number of q in and out grid points

  // variables for the fillwt function 
  int idmin, idmax, nwds;

  // variables for the setalf() function.
  double as0 = 0.118, r20 = 2.0, eps;    // as0:starting value of strong interaction constant and r20: starting value of renormalization scale

  // threshold for the heavy quarks (charm, bottom and top quarks respectively).
  double q2ch = 2.13, q2bt = 4.50, q2tp = 173.0;   

  // pdf array declaration: In this array the pdfs will be stored.
  double pdf[13];

  // type of error in the interpolation functions
  int ichk = 1;

  //----------------- QCDNUM FUNCTIONS ------------------
  QCDNUM::qcinit(lun, outfile);
  //QCDNUM::setval("qmax", 1e13);           // TODO: did not work.. set bigger value of max boundary of q^2 grid 
  QCDNUM::gxmake(x, iwt, n, nxin, nxout, iosp);     // x = xmin
  QCDNUM::gqmake(qarr, wgt, nqarr, nqin, nqout);
  QCDNUM::fillwt(pdfset_type, idmin, idmax, nwds);
  QCDNUM::setord(iord);
  QCDNUM::setalf(as0, r20);

  // grid q^2 indexes of the heavy quark thresholds
  int iqch = QCDNUM::iqfrmq(q2ch);       
  int iqbt = QCDNUM::iqfrmq(q2bt);       
  int iqtp = QCDNUM::iqfrmq(q2tp);
  // set the thresholds in the VFNS (nfix = 0)
  QCDNUM::setcbt(nfix, iqch, iqbt, iqtp); 

  // q^2 grid index of the starting scale.
  int iq0  = QCDNUM::iqfrmq(q20);     
  //evolve all pdf's.
  QCDNUM::evolfg(pdfset_type,func,pdf_flavour,iq0,eps);       //TODO: Enteder esta variable eps.    

  //cout << "     mu2       x           uv         dv        ubar        dbar         gl" << endl;
  cout << setprecision(4) << scientific;

  // loop to interpolate all pdf's at different scales of q^2
  for(int iq2 = 0; iq2 < n_q2; iq2++) 
  {
    double q2_value = q2[iq2]; 
    cout << "q2_value = " << q2_value <<endl;

    // to save the output
    ofstream myfile;
    myfile.open("/opt/qcdnum-17-01-14/output/pruebaCxx_q2_" + to_string(q2_value) +  ".csv");
    myfile << "x xuv xdv xubar xdbar xgl" << endl;

    for(int ix = 0; ix < nx-1; ix++) {
      double x_value  = x[ix];
      double uv = QCDNUM::fvalxq(pdfset_type,2,x_value,q2_value,ichk);     // the indexes are according the convention given in page 53 of the qcdnum manual
      double dv = QCDNUM::fvalxq(pdfset_type,1,x_value,q2_value,ichk);     
      double ubar = QCDNUM::fvalxq(pdfset_type,-2,x_value,q2_value,ichk);  
      double dbar = QCDNUM::fvalxq(pdfset_type,-1,x_value,q2_value,ichk);  
      double gl = QCDNUM::fvalxq(pdfset_type,0,x_value,q2_value,ichk);    
      // the above is the same as 
      // QCDNUM::allfxq(pdfset_type,x_value,q_value,pdf,0,ichk); and printing : pdf[4] << " " << pdf[5]  << " " << pdf[6]  << " " << pdf[7] << " " << pdf[8]
      
      myfile << x_value << " " << uv << " " << dv << " " << ubar << " " << dbar << " " << gl << endl;
    } 
    myfile.close();
  }                

  return 0;
}