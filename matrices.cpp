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


void free_(double **m, int n)
/* Désallocation de la mémoire pour une matrice*/
{
  int i;
  for (i=0; i<n; i++)
    {
      free(m[i]);
    }
  free(m);
}

void cholesky(double **A,int n,  double **L, double **Lt)
{
  int i,j,k,l,p;
  
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

  // Transposée de L pour obtenir Lt
  for (i=0; i<n; i++)
    {
      for (j=0; j<n; j++)
	{
	  L[i][j] = Lt[j][i];
	}
    }

}


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



  fstream f_reflexions;
  f_reflexions.open("reflexions.txt", ios::app); // à cause de ios::app, il faut supprimer à chaque fois le fichier "reflexions.txt"

  f_reflexions << "itération " << k+1 << endl; 
  f_reflexions << endl; 

  int l,p; 
  for (l=0; l<n  ; l++)
  {
    for(p=0; p<n  ; p++)
    {
      f_reflexions<< P[l][p] << " "; 
    }
    f_reflexions<< endl;
  }
  f_reflexions<< endl;
  f_reflexions<< endl; 

  f_reflexions.close(); 

  
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

    

    //fstream f4;
    //f4.open("iterations_R_TMP.txt", ios::app);
    
    //for (l=0; l<n  ; l++)
    //{
    //	for(p=0; p<n  ; p++)
    //	  {
    //	    f4<< R_TMP[l][p] << " "; 
    //	  }
    //	f4<< endl;
    //}
    //f4<< endl;
    //f4 << endl;
    //f4.close(); 

    free(x);
  }

  copie_mat(Q_TMP, Q, n);



  // Affiche les matrices Q, puis Q_t et R pour voir où est le problème
  
  //fstream f1;
  //f1.open("iterations_householder.txt", ios::app);

   
  //for (l=0; l<n  ; l++)
  //{
  // for(p=0; p<n  ; p++)
  //{
  //  f1<< Q[l][p] << " "; 
  //}
  //f1<< endl;
  //}
  //f1<< endl;

  mat_prod(Q, A, R, n, n, n);
  mat_inv(Q, TMP, n);  
  copie_mat(TMP, Q, n);


  //for (l=0; l<n  ; l++)
  //{
  //for(p=0; p<n  ; p++)
  //{
  //  f1<< Q[l][p] << " "; 
  // }
  //f1<< endl;
  //}
  //f1<< endl;


  //for (l=0; l<n  ; l++)
  //{
  //for(p=0; p<n  ; p++)
  //{
  //  f1<< R[l][p] << " "; 
  // }
  //f1<< endl;
  //}
  //f1<< endl;
  //f1<<endl;
  //f1<<endl; 

  //f1.close(); 
  

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
  f.open("iterations_schur.txt", ios::out); // Affiche la matrice A à chaque itération, qui se rapproche d'une matrice triangulaire supérieure

  f << "Au départ A = C, où C est la matrice calculée à partir de la décomposition de Cholesky" << endl;
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

    f << "Puis on affiche les matrices A à chaque itération de la décomposition de Schur" << endl;
    f << endl; 
  
  
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
