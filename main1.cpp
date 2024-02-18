#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include"bibli_fonctions.h" 
// #include"bibli_perso.h"

using namespace std;

const int N = 100; // nombre de points de la grille
const double L = 10; // Longueur de la cavité -> plus tard il faudra demander à l'utilisateur de rentrer une valeur et le programme devra s'adapter


/* Résolution de l'équation d'onde dans une cavité à 1D pour un champ électrique E(z,t) se propageant à une direction selon z. Si l'on considère une onde monochromatique l'équation devient l'équation de Helmholtz, beaucoup plus simple à résoudre. */ 




void laplacien_mat(double **D, double *z)
/* Le calcul des élements de matrice de l'opérateur laplacien en 1D peut s'effectuer à la main
RQ: Il n'y a pas besoin d'implémenter les fonctions de base (ici les fonctions tentes) tant que les 
éléments de matrice peuvent se calculer à la main */
{
  int i;
  for(i=1; i<N-1;i++)
    {
      D[i][i] = - (double(1)/(z[i+1] - z[i]) + double(1)/(z[i] - z[i-1]));
      D[i][i+1] = double(1)/(z[i+1] - z[i]);
      D[i-1][i] = double(1)/(z[i] - z[i-1]);
    }
  D[0][1] = double(1)/(z[1] - z[0]);
  D[N-2][N-1] = double(1)/(z[N-1] -z[N-2]);
  // RQ: Je ne connais pas D[0][0] et D[N-1][N-1] mais je ne m'en servirai pas
  // dans mon système matriciel!
}



void masse_mat(double **M, double *grille)
/* Le calcul des éléments de le matrice de masse peut s'effectuer à la main */
{
  M[0][0] = (grille[1] - grille[0])/3;
  M[N-1][N-1] = (grille[N-1] - grille[N-2])/3;
  int i;
  for(i=1; i<N-1; i++)
    {
      M[i][i-1] = (grille[i] - grille[i-1])/6;
      M[i][i] = (grille[i+1] - grille[i-1])/3;
      M[i][i+1] = (grille[i+1] - grille[i])/6;
    }
  M[0][1] = (grille[1] - grille[0])/6;
  M[N-1][N-2] = (grille[N-1] - grille[N-2])/6;
	  
}



void init(double *grille, double *E)
  /* La fonction init permet d'initialiser le problème: on initialise le vecteur contenant les valeurs de E selon z ainsi que la grille définissant le problème en fonction de la longueur L de la cavité; on initialise les fonctions tentes*/
{
  int i;
  double h = L/(N-1);
  for(i=0; i<N; i++)
    {
      grille[i] = h*i; // Initialisation de la grille 
    }

  E[0] = 0;
  E[N-1] = 0; // On impose les conditions aux limites de Dirichlet
  // RQ: on ne considère pas des matrices par bloc, mais on peut extraire la sous-matrice correspondant
  // à  l'action de l'opérateur sur les points intérieurs de la grille

  // Initialisation du reste des valeurs de E
  for(i=1;i<N-1; i++)
    {
      E[i] = 0;
    }

}


void alloc_(double **m, int n)
{
  int i;
  m = (double **)malloc(n*sizeof(double));
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




int main()
{
 
  double *grille = (double*)malloc(N*sizeof(double)); // subdivision de l'axe de la cavité à 1D
  // Pour l'instant la subdivision est régulière, mais on pourrait imaginer tout type de géométrie pour la grille 
  double *E = (double*)malloc(N*sizeof(double)); // vecteur contenant les valeurs du champ E 

  double** D = NULL; // élements de matrice du laplacien en 1D
  double **M = NULL; // élements de la matrice de masse en 1D 
  alloc_(M,N);
  alloc_(D,N); 


  free(grille);
  free(E);
  free_(M, N);
  free_(D, N);

  
  return 0;
}
