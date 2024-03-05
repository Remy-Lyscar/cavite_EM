#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <time.h>
#include"bibli_fonctions.h" // A decommenter lorsque je veux faire des plots Python sur les ordis du magistère (Python.h pas reconnu par VSC)
// #include"bibli_perso.h"

using namespace std;

const int N = 10; // nombre de points de la grille
const double L = 10; // Longueur de la cavité -> plus tard il faudra demander à l'utilisateur de rentrer une valeur et le programme devra s'adapter


/* Résolution de l'équation d'onde dans une cavité à 1D pour un champ électrique E(z,t) se propageant à une direction selon z. Si l'on considère une onde monochromatique l'équation devient l'équation de Helmholtz, beaucoup plus simple à résoudre. */ 




void laplacien_mat(double **D, double *z)
/* Le calcul des élements de matrice de l'opérateur laplacien en 1D peut s'effectuer à la main
RQ: Il n'y a pas besoin d'implémenter les fonctions de base (ici les fonctions tentes) tant que les 
éléments de matrice peuvent se calculer à la main */
// On ne retourne que la sous-matrice de dimension N-2 qui interviendra dans le système matriciel, sans se 
// préoccuper des éléments de matrice au niveau des bords 
{
  int i;
  for(i=0; i<N-3;i++)
    {
      D[i][i] = - (double(1)/(z[(i+1)+1] - z[i+1]) + double(1)/(z[i+1] - z[(i+1)-1])); // Décalage de 1 dans les indices
      D[i][i+1] = double(1)/(z[(i+1)+1] - z[i+1]);
      D[i+1][i] = D[i][i+1];
    }
    D[N-3][N-3] = - (double(1)/(z[N-1] - z[N-2]) + double(1)/(z[N-2] - z[N-3]));
}



void masse_mat(double **M, double *z)
/* Le calcul des éléments de le matrice de masse peut s'effectuer à la main */
{
  int i;
  for(i=0; i<N-3; i++)
    {
      M[i][i] = (z[(i+1) + 1] - z[(i+1)-1])/3;
      M[i][i+1] = (z[(i+1)+1] - z[i+1])/6;
      M[i+1][i] = (z[(i+1)+1] - z[i+1])/6;
    }	  
    M[N-3][N-3] = (z[N-1] - z[N-3])/3;
}



void init(double *grille, double *E)
  /* La fonction init permet d'initialiser le problème: on initialise le vecteur contenant les valeurs de E selon z ainsi que la grille définissant le problème en fonction de la longueur L de la cavité; 
  pas besoin d'initialiser les fonctions de base, vu que tous les calculs d'éléments de matrice peuvent s'effectuer à la main*/
{
  int i;
  double h = L/(N-1);
  for(i=0; i<N; i++)
    {
      grille[i] = h*i; // Initialisation de la grille 
    }

  // Initialisation des valeurs de E
  for(i=0;i<N-2; i++)
    {
      E[i] = 0;   // On ne crée que le vecteur des valeurs intérieurs de E, de dimension N-2
    }

}


void alloc_(double **m, int n)
{
  int i;
  m = (double **)malloc(n*sizeof(double*));
  for(i=0; i<n; i++)
    {
      m[i] = (double *)malloc(n*sizeof(double));
    }

}


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


void affichage(double *E, double* grille, double**D, double **M)
{
  fstream f; 
  f.open("verif.txt", ios::out); 

  int k; 
  for(k=0; k<N; k++)
  {
    f << grille[k] << " "; 
  }
  f << endl; 

    for(k=0; k<N-2; k++)
  {
    f << E[k] << " "; 
  }
  f << endl; 

  int i,j; 
  for (i=0; i<N-2; i++)
  {
    for(j=0; j<N-2; j++)
    {
      f << M[i][j] << " "; 
    }
    f << endl;
  }
  f << endl;

   for (i=0; i<N-2; i++)
  {
    for(j=0; j<N-2; j++)
    {
      f << D[i][j] << " "; 
    }
    f << endl;
  }
  f << endl;

  f.close(); 
}



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
  // ou de la bibliothèque du magistère

  // Pour l'instant je le fais à la main
  for (i=0; i<n; i++)
    {
      for (j=0; j<n; j++)
	{
	  L[i][j] = Lt[j][i];
	}
    }
      

  // il faudra optimiser plus tard (notamment le calcul des racines est fait deux fois
}



void aff_cholesky(double **Lt, double **L)
{
  fstream f;
  f.open("cholesky.txt", ios::out);
  
  int i,j; 
  for (i=0; i<N-2; i++)
  {
    for(j=0; j<N-2; j++)
    {
      f << Lt[i][j] << " "; 
    }
    f << endl;
  }
  f << endl;


  for (i=0; i<N-2; i++)
  {
    for(j=0; j<N-2; j++)
    {
      f << L[i][j] << " "; 
    }
    f << endl;
  }
  f << endl;

  f.close();
}



void C_calc(double **L, double **D, double **Lt, double**C, n)
{
  // On utilise une fonction de la bibliothèque du magistère

  int i;
  double **Lt_inv = (double **)malloc(n*sizeof(double*));
  for(i=0; i<n; i++)
    {
      Lt_inv[i] = (double *)malloc(n*sizeof(double));
    }

  double **L_inv = (double **)malloc(n*sizeof(double*));
  for(i=0; i<n; i++)
    {
      L_inv[i] = (double *)malloc(n*sizeof(double));
    }
  
  double **TMP = (double **)malloc(n*sizeof(double*));
  for(i=0; i<n; i++)
    {
      TMP[i] = (double *)malloc(n*sizeof(double));
    }
  mat_inv(Lt, Lt_inv, n);  // Rq: Ici on ne se sert du coup pas des propriétés des matrices L et Lt pour inverser les matrices plus rapidement ?? 

  mat_prod(D, Lt_inv, TMP, n, n, n);
  mat_prod(L_inv, TMP, C, n, n, n);

}

  

    

int main()
{
 
  double *grille = (double*)malloc(N*sizeof(double)); // subdivision de l'axe de la cavité à 1D
  // Pour l'instant la subdivision est régulière, mais on pourrait imaginer tout type de géométrie pour la grille 
  double *E = (double*)malloc((N-2)*sizeof(double)); // vecteur contenant les valeurs intérieures du champ E 

  double** D = (double **)malloc((N-2)*sizeof(double *)); // élements de matrice du laplacien en 1D
  double **M = (double **)malloc((N-2)*sizeof(double *)); // élements de la matrice de masse en 1D 


  int i;
  for (i=0; i<(N-2); i++)
    {
      D[i] = (double *)malloc((N-2)*sizeof(double));
      M[i] = (double *)malloc((N-2)*sizeof(double));
}


  // initialisation des tableaux E et grille, en prenant en compte les conditions de Dirichlet aux bords 
  init(grille, E);


  //Calcul des matrices des opérateurs Masse M et Laplacien D
  laplacien_mat(D, grille);
  masse_mat(M, grille); 

  // Affichage des vecteurs et matrice dans des fichiers pour vérification (ne marche pas sur mon ordi ?)
  // affichage(E, grille, D, M);


  // Résolution du système matriciel

  // Décomposition de Cholesky
  double** L = (double **)malloc((N-2)*sizeof(double *)); 
  double **Lt = (double **)malloc((N-2)*sizeof(double *));  


  for (i=0; i<(N-2); i++)
    {
      L[i] = (double *)malloc((N-2)*sizeof(double));
      Lt[i] = (double *)malloc((N-2)*sizeof(double));
    }

  cholesky(M, N-2, L, Lt);
  // aff_cholesky(Lt, L);

  // Calcul de la matrice C intervenant dans le problème aux valeurs propres
  // RQ: A la fin je cleanerai le code pour que toutes les matrices soient allouées à un seul endroit du code (hormis dans les fonctions qui utilisent des matrices TMP)
  double **C = (double **)malloc((N-2)*sizeof(double *));  


  for (i=0; i<(N-2); i++)
    {
      C[i] = (double *)malloc((N-2)*sizeof(double));
    }

  C_calc(L, D, Lt, C, N-2);  


  free(grille);
  free(E);
  free_(M, N-2); 
  free_(D, N-2);

  
  return 0;
}
