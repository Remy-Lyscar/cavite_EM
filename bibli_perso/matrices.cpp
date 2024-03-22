#ifndef DEF_MATRICES
#define DEF_MATRICES


#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <time.h>
#include "bibli_fonctions.h"
#include "bibli_perso.h"

using namespace std;




void cholesky(double **A,int n,  double **L, double **Lt)
{

  double **S = (double **)malloc(n*sizeof(double*)); // Sert à ne pas modifier la matrice A (on verra plus tard)
  int i,j,k,l,p;
  for (i=0; i<n; i++)
    {
      S[i] = (double*)malloc(n*sizeof(double));
    }
  
  for (i=0; i<n; i++)
    {
      Lt[i][i] = sqrt(A[i][i]);
      for (j=i+1; j<n; j++)
	{
	  Lt[i][j] = A[i][j]/sqrt(A[i][i]);
	}
      for (k=i+1; k<n; k++)
	{
	  
	  for (l=k; l<n; l++)
	    {
	      for (p=k; p<n; p++)
		{
		  A[k][l] -= Lt[i][k]*Lt[i][p];
		    }
	    }
	}
    }

  // Transposer ensuite pour avoir la matrice L (utiliser des fonctions de la bibliothèque standard
  // ou de la bibliothèque du magistère)

  // Pour l'instant je le fais à la main
  for (i=0; i<n; i++)
    {
      for (j=0; j<n; j++)
	{
	  L[i][j] = Lt[j][i];
	}
    }
      
  free_(S, n);

  // il faudra optimiser plus tard (notamment le calcul des racines est fait deux fois)
}



// à mettre dans bibli_perso
void produit_scalaire(double *a, double *b, double *rst, int n)
{
  int i; 
  *rst = 0; 
  for (i=0; i<n; i++)
  {
    *rst += a[i]*b[i];

  }
}


void norme(double *a, double *rst, int n)
{
  int i; 
  *rst = 0; 
  for (i=0; i<n; i++)
  {
    *rst += a[i] * a[i];

  }
  *rst = sqrt(*rst);
}



void vecteur_canonique(double *b, int n)
// Renvoie le vecteur b1 = (1, 0, ... ,0) de la base canonique
// Utilisé par la transformation de Householder pour la décomposition QR 

{
  int i; 
  b[0] = 1; 
  for (i=1; i<n; i++)
  {
    b[i] = 0;
  }
}


void copie_mat(double **A, double **B, int n)
// copie la matrice A dans la matrice B 
{
  int i,j; 
  for (i=0; i<n; i++)
  {
    for (j=0; j<n; j++)
    {
      B[i][j] = A[i][j];
    }
  }
}


void reflexion(double *x, double **P, int n, int k)
// Calcul de la matrice de réflexion à chaque itération de la procédure de Householder 
// Seul le bloc inférieur de dimension (n-k) est différent de l'identité

{
  int i, j, m; 
  
  double  *u= NULL; // contient le vecteur u à partir duquel on définit e
  u = (double *)malloc((n-k)*sizeof(double));

  double  *b= NULL; 
  b = (double *)malloc((n-k)*sizeof(double));

  double  *e= NULL; 
  e = (double *)malloc((n-k)*sizeof(double));

  double norme_x = 0; 
  norme(x, &norme_x, n-k); 
  vecteur_canonique(b, n-k);

  for (i=0; i<(n-k); i++)
  {
    u[i] = x[i] - norme_x*b[i]; 
  }

  double norme_u = 0; 
  norme(u, &norme_u, n-k); 

  for (i=0; i<(n-k); i++)
  {
    e[i] = u[i] / (norme_u);  
  }

  // Calcul de la matrice de réflexion


  for (m=0; m<k; m++)
    {
      for (i=0; i<n; i++)
	{
	  P[i][m] = 0;
	}
      for (j=0; j<n; j++)
	{
	  P[m][j] = 0;
	}
    }

  for (i=0; i<n; i++) // boucle inutile mais pour être sur que tous les coefs sont réinitialisés 
    {
      for (j=0; j<n; j++)
	{
	  P[i][j] = 0;
	}
    }
  
  

  for (i=0; i<k; i++)
  {
    P[i][i] = 1; 
  }

  
  for (i=0; i<(n-k); i++)
  {
    for (j=0; j<(n-k); j++)
    {
      if (j==i)
      {
        P[i+k][j+k] = 1 - 2*e[i]*e[j];
      }
      else
      {
        P[i+k][j+k] = -2*e[i]*e[j];
      }
    }
  }



  fstream f2;
  f2.open("reflexions.txt", ios::app);

  int l,p; 
  for (l=0; l<n  ; l++)
  {
    for(p=0; p<n  ; p++)
    {
      f2<< P[l][p] << " "; 
    }
    f2<< endl;
  }
  f2<< endl;
  f2<< endl; 

  f2.close(); 

  
  free(u); 
  free(e); 
  free(b); 

}


void householder(double **A, double **Q, double **R, int n)
// Effectue la décomposition QR de la matrice carrée A de taille n, en utilisant l'algorithme de Householder
{
  int i,k;
  int p,l; 

  double **Q_TMP= NULL; // Contient la matrice Q_transpose à chaque itération
  Q_TMP = (double **)malloc(n*sizeof(double *));  
  for (i=0; i<n; i++)
    {
      Q_TMP[i] = (double *)malloc(n*sizeof(double));
    }

  double **P = NULL; // Contient les matrices de réflexion à chaque itération
  P = (double **)malloc(n*sizeof(double *));  
  for (i=0; i<n; i++)
    {
      P[i] = (double *)malloc(n*sizeof(double));
    }

  double **TMP = NULL; 
  TMP = (double **)malloc(n*sizeof(double *));  
  for (i=0; i<n; i++)
    {
      TMP[i] = (double *)malloc(n*sizeof(double));
    }

  double **R_TMP = NULL; // contient la matrice R à chaque itération
  R_TMP = (double **)malloc(n*sizeof(double *));  
  for (i=0; i<n; i++)
    {
      R_TMP[i] = (double *)malloc(n*sizeof(double));
    }

  copie_mat(A, R_TMP, n); // on initialise la matrice R à A 

  
  for (k=0; k<(n-1); k++)  // PAs besoin de faire n itérations, il n'y a rien à faire sur la dernière colonne
  {
    double *x = NULL;   // contient le vecteur x à chaque itération
    x = (double *)malloc((n-k)*sizeof(double));

    for (i=k; i<n; i++)
    {
      x[i-k] = A[i][k];
    }

    reflexion(x, P, n, k);

    mat_prod(P, R_TMP, TMP, n, n, n);
    copie_mat(TMP, R_TMP, n); 

    if (k==0)
    {
      copie_mat(P, Q_TMP, n);
    }
    else
    {
      mat_prod(P, Q_TMP, TMP, n, n, n);
      copie_mat(TMP, Q_TMP, n);
    }

    fstream f4;
    f4.open("iterations_R_TMP.txt", ios::app);
    
    for (l=0; l<n  ; l++)
      {
	for(p=0; p<n  ; p++)
	  {
	    f4<< R_TMP[l][p] << " "; 
	  }
	f4<< endl;
      }
    f4<< endl;
    f4 << endl;
    f4.close(); 

    free(x);
  }

  copie_mat(Q_TMP, Q, n);

  fstream f1;
  f1.open("iterations_householder.txt", ios::app);

   
  for (l=0; l<n  ; l++)
  {
    for(p=0; p<n  ; p++)
    {
      f1<< Q[l][p] << " "; 
    }
    f1<< endl;
  }
  f1<< endl;

  mat_prod(Q, A, R, n, n, n);
  mat_inv(Q, TMP, n);  // Q: est-il vraiment nécessaire de renvoyer Q et non Q_transpose au vu de ce qu'on en fait après ?
  copie_mat(TMP, Q, n);


  for (l=0; l<n  ; l++)
  {
    for(p=0; p<n  ; p++)
    {
      f1<< Q[l][p] << " "; 
    }
    f1<< endl;
  }
  f1<< endl;


  for (l=0; l<n  ; l++)
  {
    for(p=0; p<n  ; p++)
    {
      f1<< R[l][p] << " "; 
    }
    f1<< endl;
  }
  f1<< endl;
  f1<<endl;
  f1<<endl; 

  f1.close(); 
  

  free_(Q_TMP, n);
  free_(TMP, n);
  free_(P, n);

}


void schur(double **A, double **T, double**Q, int n, int m)
// Effectue la décomposition de Schur de la matrice A, à l'aide la décomposition QR 
// réalisée par la procédure de Householder
// m est le nombre d'itérations, normalement on est censé arriver "rapidement" à une matrice 
// A_m = T triangulaire supérieure avec une précision suffisante pour en trouver les vaps et les veps 

{

  fstream f;
  f.open("iterations_schur.txt", ios::out);
  
  int i,k; 

  double **Q_TMP = NULL;  // contient les matrices Q calculées à chaque itération
  Q_TMP = (double **)malloc(n*sizeof(double *));  
  for (i=0; i<n; i++)
    {
      Q_TMP[i] = (double *)malloc(n*sizeof(double));
    } 

  double **R_TMP = NULL;  // contient les matrices R calculées à chaque itération
  R_TMP = (double **)malloc(n*sizeof(double *));  
  for (i=0; i<n; i++)
    {
      R_TMP[i] = (double *)malloc(n*sizeof(double));
    } 
  
  double **TMP = NULL; 
  TMP = (double **)malloc(n*sizeof(double *));  
  for (i=0; i<n; i++)
    {
      TMP[i] = (double *)malloc(n*sizeof(double));
    } 


  for (k=0; k<m; k++)
  {
    householder(A, Q_TMP, R_TMP, n);
    mat_prod(R_TMP, Q_TMP, TMP, n, n, n);
    copie_mat(TMP, A, n);

    
    int l,p; 
    for (l=0; l<n; l++)
      {
	for(p=0; p<n; p++)
	  {
	    f << A[l][p] << " "; 
	  }
	f << endl;
      }
    f << endl;
    f<<endl;
    

    if (m==0)
    {
      mat_inv(Q_TMP, TMP, n);   // Car on calcule en fait Q_transpose ici 
      copie_mat(TMP, Q, n);
    }
    else
    {
      mat_inv(Q_TMP, TMP, n);
      copie_mat(TMP, Q_TMP, n);
      mat_prod(Q_TMP, Q, TMP, n, n, n);
      copie_mat(TMP, Q, n);
    }
  }

  copie_mat(A, T, n); // T a priori trigsup avec suffisamment de précision pour m assez grand 

  mat_inv(Q, TMP, n);
  copie_mat(TMP,Q, n);  // Q est la matrice de passage qui contient les vecteurs propres !


  free_(Q_TMP, n);
  free_(R_TMP, n);
  free_(TMP, n);
  f.close(); 

}




#endif