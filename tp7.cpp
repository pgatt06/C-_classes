#include "tp7.hpp"
//---------------------------------------------------------------------------
//     classe Point
//---------------------------------------------------------------------------
ostream& operator<<(ostream& out , const Point& P) {
    int n= P.dim;
    if(n==0){out<<"point 0D";}
    if(n==1){out<<"("<<P.x<<")";}
    if(n==2){out<<"("<<P.x<<","<<P.y<<")";}
    if(n>2){out<<"point de dimension supérieure à 2";}
    return (out);}

//surcharge des operateurs interne à la classe

Point& Point ::operator += (const Point &Q){
    if (dim!=Q.dim){cout<<"dimension incompatibles";exit(-1);}
    if(dim==0){return (*this);}
    if (dim==1){x+=Q.x; return(*this);}
    if (dim==2){x+=Q.x; y += Q.y; return (*this);}
    else{cout<<"erreur dimension sup à 2"; exit(-1);}
};

Point& Point ::operator -= (const Point &Q){
    if (dim!=Q.dim){cout<<"dimension incompatibles";exit(-1);}
    if(dim==0){return (*this);}
    if (dim==1){x-=Q.x; return(*this);}
    if (dim==2){x-=Q.x; y -= Q.y; return (*this);}
    else{cout<<"erreur dimension sup à 2"; exit(-1);}
};

Point& Point ::operator *= (const Point &Q){
    if (dim!=Q.dim){cout<<"dimension incompatibles";exit(-1);}
    if(dim==0){return (*this);}
    if (dim==1){x+=Q.x; return(*this);}
    if (dim==2){x+=Q.x; y *= Q.y; return (*this);}
    else{cout<<"erreur dimension sup à 2"; exit(-1);}
};

Point& Point ::operator /= (const Point &Q){
    if (dim!=Q.dim){cout<<"dimension incompatibles";exit(-1);}
    if(dim==0){return (*this);}
    if (dim==1){x+=Q.x; return(*this);}
    if (dim==2){x+=Q.x; y /= Q.y; return (*this);}
    else{cout<<"erreur dimension sup à 2"; exit(-1);}
};

Point& Point::operator=(const Point& Q) {
    if (this != &Q) {
        dim = Q.dim;
        if(dim==0){return (*this);}
        if (dim==1){x=Q.x; return(*this);}
        if (dim==2){x=Q.x; y= Q.y; return (*this);}
    }
    return *this;
}


//surcharge des opérateurs (externe à la classe )
Point operator + (const Point& P, const Point & Q)
    {Point res(P); //initialisation du résultat à P 
    if(P.dim!=Q.dim){cout<<"erreur dimension";
    exit(-1);}
    else{return (res+=Q);}
    }

Point operator - (const Point& P, const Point & Q)
    {Point res(P); //initialisation du résultat à P 
    if(P.dim!=Q.dim){cout<<"erreur dimension";
    exit(-1);}
    else{return (res-=Q);}
    }
 //a definir dans les deux sens ie a*P et P*a
Point operator * (const Point& P, double a)
    {Point res(P); //initialisation du résultat à P 
    if(P.dim==0){return (res);}
    if (P.dim==1){res.x*=a; return(res);}
    if (P.dim==2){res.x*=a; res.y*=a; return (res);}
    else{cout<<"erreur dimension sup à 2"; exit(-1);}
    
    }
Point operator * (double a, const Point &P)
     {Point res(P); //initialisation du résultat à P 
    if(P.dim==0){return (res);}
    if (P.dim==1){res.x*=a; return(res);}
    if (P.dim==2){res.x*=a; res.y*=a; return (res);}
    else{cout<<"erreur dimension sup à 2"; exit(-1);}
   
    }
Point operator / (const Point& P, double a)
     {Point res(P); //initialisation du résultat à P 
     if(a==0){cout<<"division par 0";
     exit(-1);}
    if(P.dim==0){return (res);}
    if (P.dim==1){res.x/=a; return(res);}
    if (P.dim==2){res.x/=a; res.y/=a; return (res);}
    else{cout<<"erreur dimension sup à 2"; exit(-1);}
    
    }
double operator |(const Point & P, const Point & Q) //produit scalaire
    {
        if(P.dim!=Q.dim){cout<<"dimensions incompatibles"<<endl;
        exit(-1);}
        double res =0; 
        if(P.dim==0){cout<<"dimension nulle";return (0);}
        if (P.dim==1){res=P.x*Q.x;return(res);}
        if (P.dim==2){res=P.x*Q.x+P.y*Q.y; return (res);}
        else{cout<<"erreur dimension sup à 2"; exit(-1);}
    }

double norme (const Point & P){
    double res =P|P;
    return(sqrt(res));
}

Point operator * (const Point &Q , const Point &P)
     {Point res(P); //initialisation du résultat à P 
    if(P.dim==0){return (res);}
    if (P.dim==1){return (P|Q);}
    if (P.dim==2){res.x=Q.x*P.y-Q.y*P.x; res.y =-res.x; return (res);}
    else{cout<<"erreur dimension sup à 2"; exit(-1);}
    }

double mes (const Point &A, const Point &B, const Point &C)
    { double res=0;
        if(A.dim==B.dim && B.dim==C.dim){
            Point prod=((B-A)*(C-A));
            res=0.5*norme(prod);
            return(res);
        };
    exit(-1);
    }
//---------------------------------------------------------------------------
//     classe Maillage
//---------------------------------------------------------------------------
        
// Surcharge de l'opérateur << pour la classe Maillage
ostream& operator<<(ostream& out, const Maillage& m) {
    m.print(out);
    return out;
}

//---------------------------------------------------------------------------
//     classe Maillage1D
//---------------------------------------------------------------------------
Maillage1D::Maillage1D(double a, double b, int m):Maillage(1,"segments"){
    noeuds.resize(m+1);
    for(int i =0 ;i<=m;i++){
        noeuds[i]=a+i+(b-a)/(m);
    }
    numelts.resize(m);
    for(int i =0 ;i<m;i++){
        numelts.push_back({i, i + 1});
    }
};



  
//---------------------------------------------------------------------------
//     classe Maillage2D
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
//    classe EF1D
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
//    classe EF2D
//---------------------------------------------------------------------------