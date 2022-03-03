#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include "QCDNUM/QCDNUM.h"
using namespace std;

// Funciones partonicas:
//----------------------------------------------------------------------
double xupv(double x)
{
  double au = 10.19304899633023;
  double bu = 0.76, cu = 4.6, eu = 2.6, fu = 0.35, gu = 0.049;
  double pd = au * pow(x, bu) * pow((1 - x), cu) * (1 + eu * pow(x, 2) + fu * log(x) + gu * pow(log(x), 2));
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
  double A_g = -0.927244357490440, B_g = -0.52, C_g = 4.5, F_g = 0.217, G_g = 0.0112;
  double pd = A_g * pow(x, B_g) * pow((1 - x), C_g) * (1 + F_g * log(x) + G_g * pow(log(x), 2));
  return pd;
}
//----------------------------------------------------------------------
double xdbar(double x)
{
  double A_dbar = 0.14, B_dbar = -0.33, C_dbar = 24, D_dbar = 28, F_dbar = 0.071;
  double pd = A_dbar * pow(x, B_dbar) * pow((1 - x), C_dbar) * (1 - D_dbar * x + F_dbar * log(x));
  return pd;
}
//----------------------------------------------------------------------
double xubar(double x)
{
  double A_ubar = 0.14, B_ubar = -0.33, C_ubar = 11, D_ubar = 18, F_ubar = 0.071;
  double pd = A_ubar * pow(x, B_ubar) * pow((1 - x), C_ubar) * (1 - D_ubar * x + F_ubar * log(x));
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
** PDFs  xfi(x) a escala inicial iq0. iq0 declarado mas abajo.
*/
double func(int *ipdf, double *x)
{
  int i = *ipdf;
  double xb = *x;
  double f = 0;
  if (i == 0)
    f = xglu(xb); // gluon
  if (i == 1)
    f = xdnv(xb); // valence down
  if (i == 2)
    f = xupv(xb); // valence up
  if (i == 3)
    f = 0; // valence strange
  if (i == 4)
    f = xdbar(xb); // Sea quarks:
  if (i == 5)
    f = xubar(xb);
  if (i == 6)
    f = xsbar(xb);
  if (i == 7)
    f = 0; // charm de valencia?
  if (i == 8)
    f = 0; // anti-charm
  if (i == 9)
    f = 0; // bottom
  if (i == 10)
    f = 0; // anti-bottom
  if (i == 11)
    f = 0; // top
  if (i == 12)
    f = 0; // anti-top
  return f;
}

int main()
{

  // unpolized dataset, NNLO, VFNS
  int dataset_type = 1, iord = 3, nfin = 0;

  // x, q, and quarks flavour arrays
  double x[] = {1.0e-5, 2.0e-5, 3.0e-5, 5.0e-5, 6.0e-5, 7.0e-5, 8.0e-5, 9.0e-5, 
                1.0e-4, 2.0e-4, 3.0e-4, 5.0e-4, 6.0e-4, 7.0e-4, 8.0e-4, 9.0e-4, 
                1.0e-3, 2.0e-3, 3.0e-3, 5.0e-3, 6.0e-3, 7.0e-3, 8.0e-3, 9.0e-3,
                1.0e-2, 2.0e-2, 3.0e-2, 5.0e-2, 6.0e-2, 7.0e-2, 8.0e-2, 9.0e-2,
                1.0e-1, 1.5e-1, 2.0e-1, 3.0e-1, 4.0e-1, 5.5e-1, 6.0e-1, 7.0e-1, 9.0e-1};

  double q[] = {2.0e0, 2.7e0, 3.6e0, 5.0e0, 7.0e0, 1.0e1, 1.4e1,
                2.0e1, 3.0e1, 5.0e1, 7.0e1, 1.0e2, 2.0e2, 5.0e2, 1.0e3,
                3.0e3, 1.0e4, 4.0e4, 2.0e5, 1.0e6};

  // Quarks flavour composition: is an input for evolfg:
  double pdf_flavour[] =
      // tb  bb  cb  sb  ub  db   g   d   u   s   c   b   t
      {0., 0., 0., 0., 0., -1., 0., 1., 0., 0., 0., 0., 0., // 1=dval
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

  // size of x and y arrays, and min(x)
  int nx = size(x), nq = size(q);
  double qmin = q[0], qmax = q[nq - 1];    //qmin is the starting scale

  //----------------- variables for the qcinit function ------------------
  int lun = 6;
  string outfile = " "; // if lun=6, qcdnum messages appear in the standard output

  //----------  variables the gxmake function (x-grid parameters) -----------
  double xmin[] = {x[0]};
  // nxin: number of grid points. iosp: order for spline interpolation (3=quadratic)
  int iwt[] = {1}, n = 1, nxin = 100, iosp = 3; // TODO: entender iwt y n.
  int nxout;

  //----------  variables the gqmake function (q-grid parameters) -----------
  // qarr: lower and upper end of the q grid
  double qarr[] = {qmin, qmax}, wgt[] = {1e0, 1e0}; // TODO: enteder wgt
  //  number of q grid points and out grid points
  int nqin = 60, nqout;

  //----------  variables the fillwt function -----------  //TODO: entender esto de los pesos
  int id1, id2, nwds;

  //----------  variables for iqfrmq functions -----------
  double q2c = 3, q2b = 25; // thresholds, mu20  //TODO: entender esto! los thresholds

  // as0: constante de interaccion fuerte y r20: escala de renormalizacion
  double as0 = 0.364, r20 = 2.0, eps; 

  // pdf array declaration: In this array the pdfs will be stored.
  double pdf[13];

  // for select this type of error in interpolation functions
  int ichk = 1;

  //----------------- QCDNUM FUNCTIONS ------------------
  QCDNUM::qcinit(lun, outfile);
  QCDNUM::gxmake(xmin, iwt, n, nxin, nxout, iosp);
  QCDNUM::gqmake(qarr, wgt, size(qarr), nqin, nqout);
  QCDNUM::fillwt(dataset_type, id1, id2, nwds);
  QCDNUM::setord(iord);
  QCDNUM::setalf(as0, r20);

  // TODO: enteder las siguientes funciones
  int iqc = QCDNUM::iqfrmq(q2c);       // charm threshold
  int iqb = QCDNUM::iqfrmq(q2b);       // bottom threshold
  QCDNUM::setcbt(nfin, iqc, iqb, 999); // thresholds in the VFNS

  // grid index of the of the starting scale qmin (Requiered as a parameter of the evolfg function)
  int iq0  = QCDNUM::iqfrmq(qmin);     
  
  //evolve all pdf's.
  QCDNUM::evolfg(dataset_type,func,pdf_flavour,iq0,eps);       //TODO: Enteder esta variable eps.


  cout << "     mu2       x           uv         dv        ubar        dvar         gl" << endl;
  cout << setprecision(4) << scientific;

  // loop to interpolate all pdf's at different scales of q
  for(int iq = 0; iq < nq; iq = iq+1) 
  {
    double q_value = q[iq]; 
    cout << endl;

    // to save the output in a file saved in the output file
    ofstream myfile;
    myfile.open("/opt/qcdnum-17-01-14/output/pruebaCxx_q_" + to_string(q_value) +  ".csv");
    myfile << "x uv dv ubar dvar gl" << endl;

    for(int ix = 0; ix < nx; ix = ix+1) {
      double x_value  = x[ix];
      double uv = QCDNUM::fvalxq(dataset_type,2,x_value,q_value,ichk);     //TODO: aqui uso sumfxq, fvalxq o allfxq??
      double dv = QCDNUM::fvalxq(dataset_type,1,x_value,q_value,ichk);
      double ubar = QCDNUM::fvalxq(dataset_type,-2,x_value,q_value,ichk);
      double dbar = QCDNUM::fvalxq(dataset_type,-1,x_value,q_value,ichk);
      double gl = QCDNUM::fvalxq(dataset_type,0,x_value,q_value,ichk);
      cout << q_value << " " << x_value << " " << uv << " " << dv << " " << ubar << " " << dbar << " " << gl << endl;
      myfile << x_value << " " << uv << " " << dv << " " << ubar << " " << dbar << " " << gl << endl;
    } 
    myfile.close();
  } 

  cout << endl;

    // Returns all flavour-pdf values in one call. Se asigna valor a la variable pdf. Es un array de tamaño 13 de doubles  
    //QCDNUM::allfxq(dataset_type,x,q_evolution,pdf,0,1);                    

  return 0;
}