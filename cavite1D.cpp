#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <time.h>
#include "bibli_fonctions.h"
// #include"bibli_perso.h"

using namespace std;

const int N = 10; // nombre de points de la grille
const double L = 10; // Longueur de la cavité -> plus tard il faudra demander à l'utilisateur de rentrer une valeur et le programme devra s'adapter
const int MODE = 3; // mode que l'on souhaite afficher dans Python (comme comparaison des solutions analytiques et calculées)


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




// Rq: cette fonction d'allocation ne fonctionne pas (en fait je modifie le pointeur donné en argument mais pas ce 
// qu'il y a à l'adresse du pointeur) -> revoir le cours sur les pointeurs pour comprendre
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
 

void C_calc(double **L, double **D, double **Lt, double**C, int n)
// Calcul de la matrice C intervenant dans le problème aux vaps
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
  mat_inv(L, L_inv, n);

  mat_prod(D, Lt_inv, TMP, n, n, n);
  mat_prod(L_inv, TMP, C, n, n, n);

  free_(TMP, n);
  free_(L_inv, n); 
  free_(Lt_inv, n);
}






void vaps(double *V, double **T, int n)
// Renvoie les valeurs propres de T triées dans l'ordre croissant (pour l'instant j'utilise une fonction que j'ai codé qui n'est pas efficace)
{
  int i; 
  for (i=0; i<n; i++)
  {
    V[i] = T[i][i]; 
  }
}
  



int main()
{
 

  // ********** Initialisation du problème **********


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


  // ********** Résolution du système matriciel **********

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
  //aff_c(C, N-2); 
  

  // Algorithme QR appliqué à la matrice C pour déterminer ses vaps
  // RQ: la décomposition LU est plus efficace pour inverser une matrice, mais la décomposition QR 
  // est plus adaptée à la recherche de vaps et de veps

  int m = 10; // nombre d'itération de l'algorithme QR 

  double *V = (double*)malloc((N-2)*sizeof(double));  // Vecteur qui va accueillir les valeurs propres de C 

  double **T = (double **)malloc((N-2)*sizeof(double *));  
  for (i=0; i<(N-2); i++)
    {
      T[i] = (double *)malloc((N-2)*sizeof(double));
    }


  double **Q = (double **)malloc((N-2)*sizeof(double *));  
  for (i=0; i<(N-2); i++)
    {
      Q[i] = (double *)malloc((N-2)*sizeof(double));
    }

  schur(C, T, Q, N-2, m);
  aff_T(T, N-2);
  vaps(V, T, N-2);

  // On en déduit les vecteurs d'onde k possibles

  double *k = (double*)malloc((N-2)*sizeof(double));
  for (i=0; i<(N-2); i++)
  {
    k[i] = sqrt(-(V[i]));
  }
  aff_k(k, N-2);
  // Pour chaque valeur propre k_i le vecteur propre associé est la i-eme colonne de la matrice Q de la décomposition de Schur

  // D'abord on calcul uniquement le vep associé au mode 3 pour avoir un beau plot, ensuite on fera un truc plus élégant
  
  for (i=0; i<N-2; i++)
    {
      E[i] = Q[i][5];
    }

  fstream vep;
  vep.open("vep.dat", ios::out);
  for (i=0; i<N-2; i++)
    {
      vep << E[i] << endl;
    }
  vep.close();


  fstream fgrille;
  fgrille.open("grille.dat", ios::out); 
  for(i=0; i<N-2; i++)
  {
   fgrille << grille[i] << endl; 
  }
  fgrille.close(); 

  // ********** Affichage des résultats - Comparaison avec les solutions analytiques **********


  // Affichage de la grille et des fonctions tentes


  // Solutions analytiques (pour le mode choisi par l'utilisateur): on ne représente pour l'instant que l'enveloppe
  
  // il faut créer un fichier contenant les données du problème pour que Python puisse s'y retrouver 
  
  ostringstream pyth; 
    pyth 
    << "import numpy as np\n"
    << "Z = np.linspace(0,10, 1000)\n"
    << "Y = [np.sin((3*np.pi*z)/10) for z in Z]\n"
    << "E = loadtxt('vep.dat')\n"
    << "G = loadtxt('grille.dat')\n"
    << "for i in range(len(E)): \n"
    << "  plot(G[i], E[i], '+')\n"
    << "plot(Z,Y)\n"
    ; 

  make_plot_py(pyth);

  



  // Solutions calculées (pour le mode choisi par l'utilisateur)








  free(grille);
  free(E);
  free_(M, N-2); 
  free_(D, N-2);
  free_(L, N-2); 
  free_(Lt, N-2);
  free(V); 
  free(k); 
  free_(Q, N-2); 
  free_(T, N-2); 
  free_(C, N-2);

  
  return 0;
}
