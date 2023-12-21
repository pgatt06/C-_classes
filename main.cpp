#include "tp7.hpp"
#include <ctime>
#include <cstdlib>

double k=1;         //variable globale du main
double pi=4*atan(1);//variable globale du main
// solution et donnée Helmholtz sur [0,1]
double uex1(const Point& P){return cos(pi*P.x);}
double f1(const Point& P)  {return (k*k-pi*pi)*uex1(P);}
// solution et donnée Helmholtz sur [0,1]x[0,1]
double uex2(const Point& P){return cos(pi*P.x)*cos(pi*P.y);}
double f2(const Point& P)  {return (k*k-2*pi*pi)*uex2(P);}

int main(int argc, char * argv[] )
{
  //srand(time(0));       //initialisation du générateur aleatoire
  cout<<"------------ test Point ----------------"<<endl;
  Point A(4), B(1), C(3);
  cout<<"2*(B-C)+A="<<2*(B-C)+A<<endl;
  Point O(0,0), P(0,1), Q(1,1), R(1,2);
  cout<<"O+P+2.*Q-R/2.="<<(O+P+2.*Q-R/2.) <<endl;
  cout<<"------------ test Maillage 1D ----------------"<<endl;
  Maillage1D mails(0,1,10);
  cout<<mails<<endl;
  //cout<<"------------ test Maillage2D ----------------"<<endl;
  //Maillage2D mailc(5,3);
  //cout<<mailc<<endl;
  //Maillage2D mailr(2,4,10,20,4,5);
  // cout<<mailr<<endl;
  // // resoudre (delta+k^2)u = f sur Omega, du/dn=0 sur dOmega
  // cout<<"------------ test Helmholtz 1D ----------------"<<endl;
  // Maillage1D mail1(0,1,20);
  // Helmholtz<Maillage1D,EF1D> helm1(mail1,k,f1,uex1);
  // helm1.assembleMatrices();
  // helm1.resoudre();
  // cout<<"solution Helmholtz 1D : "<<helm1.sol<<endl;
  // // cout<<"erreur L2 = "<< helm1.calculErreur()<<endl;
  // // helm1.exporte("sol1D");
  // cout<<"------------ test Helmholtz 2D ----------------"<<endl;
  // Maillage2D mail2(20,20);
  // Helmholtz<Maillage2D,EF2D> helm2(mail2,k,f2,uex2);
  // helm2.assembleMatrices();
  // helm2.resoudre();
  // // cout<<"erreur L2 = "<< helm2.calculErreur()<<endl;
  // // helm2.exporte("sol2D");
  return (0);
}
