#ifndef SPARSE_HPP
#define SPARSE_HPP
#include <iostream>
#include <map>
#include <vector>
#include <utility>
#include <cmath>
using namespace std;

//---------------------------------------------------------------------------
//    opération algébriques sur vector<T>
//---------------------------------------------------------------------------
template<typename T> vector<T> operator+(const vector<T>& u, const vector<T>& v)
{
    vector<T> w(u);
    auto itv = v.begin();
    for(auto itw=w.begin(); itw!=w.end(); ++itw, ++itv) *itw+=*itv;
    return w;
}
template<typename T> vector<T> operator-(const vector<T>& u, const vector<T>& v)
{
    vector<T> w(u);
    auto itv = v.begin();
    for(auto itw=w.begin(); itw!=w.end(); ++itw, ++itv) *itw-=*itv;
    return w;
}
template<typename T> vector<T> operator*(const vector<T>& u, const T& s)
{
    vector<T> w(u);
    for(auto& wi : w) wi*=s;
    return w;
}
template<typename T> vector<T> operator*(const T& s, const vector<T>& u)
{
    vector<T> w(u);
    for(auto& wi : w) wi*=s;
    return w;
}
template<typename T> vector<T> operator/(const vector<T>& u, const T& s)
{
    vector<T> w(u);
    for(auto& wi : w) wi/=s;
    return w;
}
template<typename T> T operator|(const vector<T>& u, const vector<T>& v)
{
    T s=0.;
    auto itv = v.begin();
    for(auto itu=u.begin(); itu!=u.end(); ++itu, ++itv) s+=*itu * *itv;
    return s;
}
template<typename T> T norme(const vector<T>&u)
{
    return sqrt(u|u);
}
template<typename T> ostream& operator<<(ostream& os,const vector<T>& v)
{
  os<<"(";
  auto itv=v.begin();
  for(;itv!=v.end()-1;++itv) os<<(*itv)<<",";
  os<<(*itv)<<")";
  return os;
}
//---------------------------------------------------------------------------
//     opérations sur Pint (pair<int,int>)
//---------------------------------------------------------------------------
typedef pair<int,int> Pint;
inline bool operator<(const Pint& ij1, const Pint& ij2)
{
  if(ij1.first<ij2.first) return true;
  if(ij1.first>ij2.first) return false;
  if(ij1.second<ij2.second) return true;
  return false;
}
inline ostream& operator<<(ostream& os,const Pint& ij)
{ os<<"("<<ij.first<<","<<ij.second<<")"; return os;}

//---------------------------------------------------------------------------
//     classe Sparse
//---------------------------------------------------------------------------
template <typename T>
class SparseT : public map<Pint,T>
{public:
  int m,n; //dimensions de la matrice
  SparseT(int mi=0,int ni=0): m(mi),n(ni){}
  T& operator()(int i, int j);
  T& operator()(int i, int j) const;
  void supprime(int i,int j);
  SparseT<T>& operator+=(const SparseT<T>& M);  // +=M
  SparseT<T>& operator-=(const SparseT<T>& M);  // -=M
  SparseT<T>& operator*=(const T& s);           // *=s
  SparseT<T>& operator/=(const T& s);           // /=s
  void remplissage() const; // affiche le remplissage
};
//fonctions associées à la classe SparseT<T>
//acces element (i,j) en lecture/ecriture
template <typename T> T& SparseT<T>::operator()(int i, int j)
{
 if(i<1 || j<1) {cout<<"indices négatifs ou nuls\n";exit(-1);}
 Pint ij(i,j);
 auto itm=this->find(ij);
 if(itm==this->end())
 {
   if(i>m) m=i;
   if(j>n) n=j;
   return (*this)[ij]=T();
 }
 return itm->second;
}
//accès element (i,j) en lecture seulement
template <typename T> T& SparseT<T>::operator()(int i, int j) const
{
 if(i<1 || j<1) {cout<<"indices négatifs ou nuls\n";exit(-1);}
 auto itm=this->find(Pint(i,j));
 if(itm==this->end()) return T();
 return itm->second;
}
//supprime l'élément (i,j) si il existe
template <typename T>void SparseT<T>::supprime(int i,int j)
{
 auto itm=this->find(Pint(i,j));
 if(itm!=this->end()) this->erase(itm);
}
// ajoute une matrice sparse à la matrice courante C+=M
template <typename T> SparseT<T>& SparseT<T>::operator+=(const SparseT<T>& M)
{
  for(auto itm=M.begin();itm!=M.end();++itm) //parcourt matrice M
  {
    auto itmc=this->find(itm->first);
    if(itmc!=this->end())
    {
      itmc->second+=itm->second;
      if(itmc->second==T()) this->erase(itmc);
    }
    else
    {
      (*this)[itm->first]=-itm->second;
      if(itm->first.first>m) m=itm->first.first;
      if(itm->first.second>m) n=itm->first.second;
    }
  }
 return *this;
}
// retire une matrice sparse à la matrice courante C-=M
template <typename T> SparseT<T>& SparseT<T>::operator-=(const SparseT<T>& M)
{
  for(auto itm=M.begin();itm!=M.end();++itm) //parcourt matrice M
  {
    auto itmc=this->find(itm->first);
    if(itmc!=this->end())
    {
     itmc->second-=itm->second;
     if(itmc->second==T()) this->erase(itmc);
    }
    else
    {
      (*this)[itm->first]=-itm->second;
      if(itm->first.first>m) m=itm->first.first;
      if(itm->first.second>m) n=itm->first.second;
    }
  }
 return *this;
}
// produit par un scalaire C*=s
template <typename T>
SparseT<T>& SparseT<T>::operator*=(const T& s)
{
 for(auto itm=this->begin();itm!=this->end();++itm) //parcourt matrice
    itm->second*=s;
 return *this;
}
// division par un scalaire C/=s
template <typename T> SparseT<T>& SparseT<T>::operator/=(const T& s)
{
 if(s==T()){cout<<"division par 0\n";exit(-1);}
 for(auto itm=this->begin();itm!=this->end();++itm) //parcourt matrice
    itm->second/=s;
 return *this;
}
//affichage du remplissage
template <typename T> void SparseT<T>::remplissage() const
{
 cout<<" ";
 for(int j=1;j<=n;j++) cout<<j%10;
 cout<<endl;
 for(int i=1;i<=m;i++)
 {
   cout<<i%10;
   for(int j=1;j<=n;j++)
   if(this->find(Pint(i,j))==this->end()) cout<<" ";
   else cout<<"x";
   cout<<endl;
 }
}
// neutre/opposé
template <typename T> SparseT<T> operator+(const SparseT<T>& M)
{ return M;}
template <typename T> SparseT<T> operator-(const SparseT<T>& M)
{ return -1.*M;}
// somme de deux matrices sparse
template <typename T> SparseT<T> operator+(const SparseT<T>& M1, const SparseT<T>& M2)
{ return SparseT<T>(M1)+=M2;}
// différence de deux matrices sparse
template <typename T> SparseT<T> operator-(const SparseT<T>& M1, const SparseT<T>& M2)
{ return SparseT<T>(M1)-=M2;}
// produit d'une matrice sparse par un scalaire
template <typename T> SparseT<T> operator*(const SparseT<T>& M, const T& s)
{ return SparseT<T>(M)*=s;}
template <typename T> SparseT<T> operator*(const T& s,const SparseT<T>& M)
{ return SparseT<T>(M)*=s;}
// division d'une matrice sparse par un scalaire
template <typename T> SparseT<T> operator/(const SparseT<T>& M, const T& s)
{ return SparseT<T>(M)/=s;}
//affichage de la matrice sparse
template<typename T> ostream& operator<<(ostream& os,const SparseT<T>& M)
{
  os<<"matrice sparse "<<M.m<<"x"<<M.n<<", nb coeff="<<M.size()<<endl;
  for(auto itm=M.begin();itm!=M.end();++itm)
      os<<"  "<<itm->first<<" : "<<itm->second<<endl;
  return os;
}
//produit matrice sparse x vector
template<typename T> vector<T> operator*(const SparseT<T>& M, const vector<T>& v)
{
  vector<T> r(M.m,T());
  for(auto itm=M.begin();itm!=M.end();++itm)
    r[itm->first.first-1]+=itm->second*v[itm->first.second-1];
  return r;
}

//resolution de Ax=b par gradient conjugué (matrice symétrique def positive)
template<typename T> vector<T> gradConj(const SparseT<T>& A, const vector<T>& b, const vector<T>& x0, int kmax, const T& eps=1e-4)
{
    vector<T> x=x0, g=A*x-b, w=g, aw;
    T nb=norme(b), ng=norme(g), ngp1, tk;
    int k=0;
    for(;k<kmax;k++)
    {
        aw=A*w;
        tk=(g|w)/(aw|w);
        x=x-tk*w;
        g=g-tk*aw;
        ngp1=norme(g);
        if(ngp1<eps*nb)  break;//convergence
        w=g+(ngp1*ngp1)/(ng*ng)*w;
        ng=ngp1;
    }
    if(k==kmax) cout<<"non convergence du gradient conjugue en "<<kmax;
    else        cout<<"convergence du gradient conjugue en "<<k;
    cout<<" iterations, residu = "<<ngp1<<endl;
    return x;
}
#endif
