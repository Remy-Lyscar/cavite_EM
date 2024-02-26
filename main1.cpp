#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <time.h>
// #include"bibli_fonctions.h" 
// #include"bibli_perso.h"

using namespace std;

const int N = 100; // nombre de points de la grille
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

    for(k=0; k<N; k++)
  {
    f << E[k] << " "; 
  }
  f << endl; 

  int i,j; 
  for (i=0; j<N-2; j++)
  {
    for(j=0; j<N-2; j++)
    {
      f << M[i][j] << " "; 
    }
    f << endl;
  }
  f << endl;

   for (i=0; j<N-2; j++)
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



int main()
{
 
  double *grille = (double*)malloc(N*sizeof(double)); // subdivision de l'axe de la cavité à 1D
  // Pour l'instant la subdivision est régulière, mais on pourrait imaginer tout type de géométrie pour la grille 
  double *E = (double*)malloc((N-2)*sizeof(double)); // vecteur contenant les valeurs intérieures du champ E 

  double** D = NULL; // élements de matrice du laplacien en 1D
  double **M = NULL; // élements de la matrice de masse en 1D 
  alloc_(M,N);
  alloc_(D,N); 


  // initialisation des tableaux E et grille, en prenant en compte les conditions de Dirichlet aux bords 
  init(grille, E);


  //Calcul des matrices des opérateurs Masse M et Laplacien D
  laplacien_mat(D, grille);
  masse_mat(M, grille); 

  // Affichage des vecteurs et matrice dans des fichiers pour vérification (ne marche pas sur mon ordi ?)
  // affichage(E, grille, D, M)



  // Résolution du système matriciel 


  free(grille);
  free(E);
  free_(M, N-2); 
  free_(D, N-2);

  
  return 0;
}
