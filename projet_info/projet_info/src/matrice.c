//
//  matrice.c
//  projet_info
//
//  Created by Alban Laborde-Laulhé on 23/09/2025.
//


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "matrice.h"


//M[i] = Pointeur de la ligne
Matrice createMatrice(int n, int m) //LiCo //Crée une matrice de taille m*n sans valeur particulière
{
    Matrice matrice;
    matrice.n = n; //ligne
    matrice.m = m; //colonne
    //Allocation du tableau de pointeur
    matrice.data = malloc(n * sizeof(double));
    
    //Allocation des lignes
    for( int i = 0; i<n ; i++)
    {
        matrice.data[i] = malloc(m * sizeof(double));
        for(int j=0; j<m; j++)
        {
            matrice.data[i][j] = 0;
        }
    }
    
    return matrice;
}

Matrice copieMatrice(Matrice M)
{
    int n = M.n, m = M.m;
    Matrice matrice;
    matrice.n = n; //ligne
    matrice.m = m; //colonne
    //Allocation du tableau de pointeur
    matrice.data = malloc(n * sizeof(double));
    
    //Allocation des lignes
    for( int i = 0; i<n ; i++)
    {
        matrice.data[i] = malloc(m * sizeof(double));
        for(int j=0; j<m; j++)
        {
            matrice.data[i][j] = M.data[i][j];
        }
    }
    
    return matrice;
    
}

void afficheMatrice(Matrice mat) //Affiche la matrice
{
    for(int i = 0; i<mat.n; i++)
    {
        for(int j=0; j<mat.m;j++)
        {
            printf("%8.2f", mat.data[i][j]);
        }
        printf("\n");
    }
}

//Génére une matrice de test où chaque coeff = id_coeff
void testMatrice(Matrice M)
{
    int k = 0;
    for(int i=0; i<M.n; i++)
    {
        for(int j=0; j<M.m; j++)
        {
            M.data[i][j] = k;
            k++;
        }
    }
}

//Crée la matrice transposé
Matrice transpose(Matrice M)
{
    double a = M.n, b = M.m;
    //Dans M : 0<i<a ; 0<j<b
    //Dans MT : 0<i<b; 0<j<a
    Matrice MT = createMatrice(b, a);
    for(int i=0; i<a; i++)
    {
        for(int j=0; j<b; j++)
        {
            MT.data[j][i] = M.data[i][j];
        }
    }
    
    return MT;
}

//Crée le produit de deux matrices
Matrice produit(Matrice M, Matrice N)
{
    if(M.m != N.n)
    {
        printf("Produit Matriciel invalide !! \n");
        exit(EXIT_FAILURE);
    }
    int n_somme = M.m;
    Matrice res = createMatrice(M.n, N.m);
    for(int i=0; i<res.n; i++)
    {
        for(int j=0; j<res.m; j++)
        {
            double coeff = 0;
            for(int k = 0; k<n_somme; k++)
            {
                coeff += M.data[i][k]*N.data[k][j];
            }
            res.data[i][j] = coeff;
        }
    }
    return res;
}

void pScaMat(Matrice M, double a){
    for(int i=0; i<M.n; i++)
    {
        for(int j=0; j<M.m; j++)
        {
            M.data[i][j] = a*M.data[i][j];
        }
    }
}

double pScalaire(Matrice X)
{
    if(X.m != 1)
    {
        printf("Pas un vecteur colonne !!! \n");
        exit(EXIT_FAILURE);
    }
    Matrice TX = transpose(X);
    return produit(X, TX).data[0][0];
    
}
//Addition M + aN
Matrice addMat(Matrice M, Matrice N, double a)
{
    if(M.n != N.n || M.m != N.m)
    {
        printf("Addition entre deux matrices de tailles différentes !! \n");
        exit(EXIT_FAILURE);
    }
    int n = M.n, m = M.m;
    Matrice res = createMatrice(n, m);
    
    for(int i=0; i<n;i++)
    {
        for(int j=0; j<m; j++)
        {
            res.data[i][j] = M.data[i][j] + a*N.data[i][j];
        }
    }
    return res;
}

double norme(Matrice M)
{
    return sqrt(pScalaire(M));
}

//Génère une matrice identité de taille n
Matrice Identite(int n)
{
    Matrice I = createMatrice(n, n);
    for(int i=0;i<n;i++)
    {
        for(int j=0; j<n; j++)
        {
            I.data[i][j] = i == j ? 1 : 0;
        }
    }
    
    return I;
}

//Trouver une valeur propre par la méthode de Jacobi

void jacobi(double **a, double **b, int n)
{
  int i,j,k,l;
  double M,alpha,t,h,c,s,r,u,v;
  int iter=0;
    
  for (i=1;i<=n;i++)
       { for (j=1;j<=n;j++) b[i][j]=0;
         b[i][i]=1;
       }

  do {
    M=0;
      for (i=1;i<n;i++)
      for (j=i+1;j<=n;j++)
        if (fabs(a[i][j])>M) {
          k=i;
          l=j;
          M=fabs(a[i][j]);
        }
         
      if (M<TOL) break;
      iter++;
        alpha=0.5*(a[l][l]-a[k][k])/a[k][l];
      t=sqrt(1+alpha*alpha);
      if (alpha<0)
        t=-t;
      t-=alpha;
      h=t*a[k][l];
      c=1/sqrt(1+t*t);
      s=t*c;
      r=s/(1+c);
        for (j=1;j<=n;j++) {
        u=b[j][k];
              v=b[j][l];
              b[j][k]-=s*(v+r*u);
              b[j][l]+=s*(u-r*v);
              if (j!=k && j!=l) {
          u=a[k][j];
                  v=a[l][j];
                  a[k][j]-=s*(v+r*u);
                  a[l][j]+=s*(u-r*v);
                  a[j][k]=a[k][j];
                  a[j][l]=a[l][j];
        }
      }
        a[k][k]-=h;
        a[l][l]+=h;
        a[k][l]=a[l][k]=0;
  } while (M>TOL && iter<1000);
  printf("iter = %d\n",iter);
}

/*
//Trouver une valeur propre par méthode de la puissance inverse
//Mu = 1 selon l'article
double vpPI(Matrice M, int maxITER)
{
    if(M.m != M.n)
    {
        printf("Recherche de valeur propre sur une matrice non carrée !! \n");
        exit(EXIT_FAILURE);
    }
    int N = M.n;
    
    Matrice b = createMatrice(N, 1);
    //On sait que les coeffs de b sont proche de 1 (article)
    for(int i=0; i<N; i++){
        b.data[i][0] = 1;
    }
    
    //Calculons la constante c0
    //Soit la Matrice M - muI = M - I
    Matrice calc = addMat(M, Identite(N), -1);
    produit(calc, b);
    double c = norme(calc);
    
    //C'est parti pour converger (linéairement) wouhouu
    for(int k=0; k<maxITER;k++)
    {
        b = pScaMat(produit(calc, b), 1/c);
    }
    
    return

}
*/



