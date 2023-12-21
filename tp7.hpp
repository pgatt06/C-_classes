#ifndef TP7_HPP
#define TP7_HPP
#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <utility>
#include "sparse.hpp"
using namespace std;
//---------------------------------------------------------------------------
//     typedef
//---------------------------------------------------------------------------
typedef vector<int> Numeros;
typedef vector<double> Vecteur;
typedef SparseT<double> Sparse;
typedef vector<vector<double>> Matrix;
//---------------------------------------------------------------------------
//     classe Point  (Point 2D ou 1D)
//---------------------------------------------------------------------------
class Point
{public:
  int dim;     //dimension du point (1 ou 2)
  double x, y; //abcisse, ordonnée
  Point(){dim =0;};                 //constructeur 0D
  Point(double a){dim=1; x=a;};          //constructeur 1D
  Point(double a, double b){dim =2; x=a; y=b;}; //constructeur 2D
  Point(const Point & P){dim=P.dim; //constructeur par copie 
  if(dim==1){x=P.x;};
  if(dim==2){x=P.x;y=P.y;}}

  //surcharge des opérateurs internes à la classe 
    Point & operator += (const Point &Q);
    Point & operator -= (const Point &Q);
    Point & operator *= (const Point &Q);
    Point & operator /= (const Point &Q);
     Point & operator = (const Point &Q);
};
ostream & operator<<(ostream & out,const Point& P);


//surcharge des opérateurs (externe à la classe )
Point operator + (const Point& P, const Point & Q);
Point operator - (const Point& P, const Point & Q);
// //a definir dans les deux sens ie a*P et P*a
Point operator * (const Point& P, double a);
Point operator * (double a, const Point &P);
Point operator / (const Point& P, double a);
double operator |(const Point & P, const Point & Q); //produit scalaire entre P et Q 
double norme (const Point & P);

typedef double (*fun)(const Point&);  // pointeur de fonction Point->double

double mes (const Point &A, const Point &B, const Point &C);
//---------------------------------------------------------------------------
//     classe Maillage
//---------------------------------------------------------------------------

class Maillage
{public:
 int dim;                                  // dimension de la géométrie (1,2)
 string geometrie;                         // nom de la géométrie maillée
 vector<Point> noeuds;                     // liste des noeuds comptés une seule fois
 list<Numeros> numelts;                    // liste des numéros des noeuds


 Maillage(int d=0, const string& geo="") : dim(d), geometrie(geo) {}


 void print(ostream& out=cout) const{
        cout<<"Maillage du"<<(dim == 1 ? "segment" : "triangle")<<geometrie<<endl;
        cout <<"liste des noeuds ( " << noeuds.size() <<"points )"<<endl;
        for (const Point& p : noeuds) {
        cout << "(" << p.x << "," << p.y << "),";
    }
        cout << endl;

       cout << "Liste des elements (" << numelts.size()/2 << " " << (dim == 1 ? "segments" : "triangles") << ")" << endl;
        for (const Numeros & N : numelts) {
            out << "(" << N << "),";
        }
        cout << endl;
    
    };;      // affichage du maillage
};

ostream& operator<<(ostream& out,const Maillage& M);
//---------------------------------------------------------------------------
//     classe Maillage1D
//---------------------------------------------------------------------------
class Maillage1D : public Maillage
{public:
  Maillage1D(double a, double b, int m);   //maillage du segment [a,b] avec m>0 intervalles;
};
//---------------------------------------------------------------------------
//     classe Maillage2D
//---------------------------------------------------------------------------
class Maillage2D : public Maillage
{public :
  Maillage2D(int m,int n); //constructeur maillage du carre [0,1]x[0,1]
  Maillage2D(double a,double b,double c,double d,int m,int n); //maillage du rectangle [a,b]x[c,d]
  void maille_carre_unite(int m,int n);                        //maillage du carré unité
};
//---------------------------------------------------------------------------
//     classe EF
//---------------------------------------------------------------------------
class EF   // classe parent
{public :
  EF(int d): dim(d){}
  virtual ~EF(){}
  int dim;  //dimension de l'élément (1 ou 2)
  virtual void masseP1(const vector<Point>&, Matrix&) const =0; // matrice de masse élémentaire
  virtual void rigidP1(const vector<Point>&, Matrix&) const =0; // matrice de rigidité élémentaire
};
//---------------------------------------------------------------------------
//     classe EF1D
//---------------------------------------------------------------------------
//...
//---------------------------------------------------------------------------
//     classe EF2D
//---------------------------------------------------------------------------
//...
//---------------------------------------------------------------------------
//     classe Helmholtz<MT,EFT>
//---------------------------------------------------------------------------
template <typename MT, typename EFT>
class Helmholtz
{protected:
  const MT* mail_ = nullptr;         //pointeur sur un maillage de type MT
  EFT ef_;                           //Elément fini de type EFT
  double k_;                         //nombre d'onde du problème
  fun f_, uex_;                      //fonction f et uex si connue
 public:
  Sparse M,K;                        //matrices K et M
  Vecteur sol;                       //vecteur solution
  Helmholtz(const MT& m, double k, fun f, fun uex=nullptr); //constructeur
  void assembleMatrices();           //calcul assemblé des matrices K et M
  Vecteur& resoudre();               //résolution du systeme (K-k^2*M)*U=M*F
  // double calculErreur() const;
  // void exporte(const string& fn) const;
};
// implementation des fonctions de la classe template Helmholtz
//...
#endif
