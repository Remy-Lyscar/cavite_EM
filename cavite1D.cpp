#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <time.h>
#include"bibli_perso.h" // contient aussi bibli_fonctions, la bibiothèque du magistère

using namespace std;

const int N = 10; // nombre de points de la grille
const double l = 30; // Longueur de la cavité 
const int MODE = 3; // Mode de la cavité que l'on veut afficher à la fin pour comparer avec la solution analytique 


/* Résolution de l'équation d'onde dans une cavité à 1D pour un champ électrique E(z,t) se propageant à une direction selon z. Si l'on considère une onde monochromatique l'équation devient l'équation de Helmholtz, beaucoup plus simple à résoudre. */ 


void laplacien_mat(double **D, double *z, int n)
/* Le calcul des élements de matrice de l'opérateur laplacien en 1D peut s'effectuer à la main
RQ: Il n'y a pas besoin d'implémenter les fonctions de base (ici les fonctions tentes) tant que les 
éléments de matrice peuvent se calculer à la main */
// On ne retourne que la sous-matrice de dimension N-2 qui interviendra dans le système matriciel, sans se 
// préoccuper des éléments de matrice au niveau des bords 
{
  int i;
  for(i=0; i<n-3;i++)
    {
      D[i][i] = - (double(1)/(z[(i+1)+1] - z[i+1]) + double(1)/(z[i+1] - z[(i+1)-1])); // Décalage de 1 dans les indices
      D[i][i+1] = double(1)/(z[(i+1)+1] - z[i+1]);
      D[i+1][i] = D[i][i+1];
    }
    D[n-3][n-3] = - (double(1)/(z[n-1] - z[n-2]) + double(1)/(z[n-2] - z[n-3]));
}



void masse_mat(double **M, double *z, int n)
/* Le calcul des éléments de le matrice de masse peut s'effectuer à la main */
{
  int i;
  for(i=0; i<n-3; i++)
    {
      M[i][i] = (z[(i+1) + 1] - z[(i+1)-1])/3;
      M[i][i+1] = (z[(i+1)+1] - z[i+1])/6;
      M[i+1][i] = (z[(i+1)+1] - z[i+1])/6;
    }	  
    M[n-3][n-3] = (z[n-1] - z[n-3])/3;
}



void init(double *grille, double *E, int n, double l)
  /* La fonction init permet d'initialiser le problème: on initialise le vecteur contenant les valeurs de E selon z ainsi que la grille définissant le problème en fonction de la longueur L de la cavité; 
  pas besoin d'initialiser les fonctions de base, vu que tous les calculs d'éléments de matrice peuvent s'effectuer à la main*/
{
  int i;
  double h = l/(n-1);
  for(i=0; i<n; i++)
    {
      grille[i] = h*i; // Initialisation de la grille 
    }

  // Initialisation des valeurs de E
  for(i=0;i<(n-2); i++)
    {
      E[i] = 0;   // On ne crée que le vecteur des valeurs intérieurs de E, de dimension N-2
    }

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
  mat_inv(Lt, Lt_inv, n);
  mat_inv(L, L_inv, n);

  mat_prod(D, Lt_inv, TMP, n, n, n);
  mat_prod(L_inv, TMP, C, n, n, n);

  free_(TMP, n);
  free_(L_inv, n); 
  free_(Lt_inv, n);
}




void vaps(double *V, double **T, int n)
// Renvoie les valeurs propres de T, et les conserve dans le vecteur V 
{
  int i; 
  for (i=0; i<n; i++)
  {
    V[i] = T[i][i]; 
  }
}

/* Classement des elements d'un tableau par ordre croissant */


void tri_croissant(double *tabl, int n)
{

  int i, j; 
  double piv; 

  for (i = 1; i < n; i++){
    piv = tabl[i];
    for (j = i-1; j >= 0; j--)
      if (piv >= tabl[j])
	break;
      else
	tabl[j+1] = tabl[j];
    tabl[j+1] = piv;
  }
}
  



int main()
{
 

  // ********** Initialisation du problème **********


  double *grille = (double*)malloc(N*sizeof(double)); // subdivision de l'axe de la cavité à 1D
  double *E = (double*)malloc((N-2)*sizeof(double)); // vecteur contenant les valeurs intérieures du champ E 
  double** D = (double **)malloc((N-2)*sizeof(double *)); // élements de matrice du laplacien en 1D
  double **M = (double **)malloc((N-2)*sizeof(double *)); // élements de la matrice de masse en 1D
  double **C = (double **)malloc((N-2)*sizeof(double *));  // matrice intervenant dans le problème aux vaps
  double** L = (double **)malloc((N-2)*sizeof(double *));  // matrices L et Lt: résultat de la décomposition de Cholesky de la matrice de masse M
  double **Lt = (double **)malloc((N-2)*sizeof(double *));
  double *V = (double*)malloc((N-2)*sizeof(double));  // Vecteur qui va accueillir les valeurs propres de C 
  double **T = (double **)malloc((N-2)*sizeof(double *));  // Matrice de Schur, contient les vaps de C à la fin du programme 
  double **Q = (double **)malloc((N-2)*sizeof(double *));  // Matrice de passage, ie contient les veps de la matrice C 


  int i;
  for (i=0; i<(N-2); i++)
    {
      D[i] = (double *)malloc((N-2)*sizeof(double));
      M[i] = (double *)malloc((N-2)*sizeof(double));
      C[i] = (double *)malloc((N-2)*sizeof(double));
      L[i] = (double *)malloc((N-2)*sizeof(double));
      Lt[i] = (double *)malloc((N-2)*sizeof(double));
      T[i] = (double *)malloc((N-2)*sizeof(double));
      Q[i] = (double *)malloc((N-2)*sizeof(double));
    }


  

  // initialisation des tableaux E et grille, en prenant en compte les conditions de Dirichlet aux bords 
  init(grille, E, N, l);


  //Calcul des matrices des opérateurs Masse M et Laplacien D
  laplacien_mat(D, grille,N);
  masse_mat(M, grille,N); 

  // Affichage des vecteurs initialisés et des matrices utilisés pour la méthode des éléments finis
  affichage(E, grille, D, M,N);


  


  // ********** Résolution du système matriciel **********

  
  // Décomposition de Cholesky de la matrice de masse
  cholesky(M, N-2, L, Lt);
  // aff_cholesky(Lt, L, N);
  
  // Calcul de la matrice C intervenant dans le problème aux valeurs propres
  C_calc(L, D, Lt, C, N-2);
  //aff_c(C, N-2); 
  

  // Algorithme QR appliqué à la matrice C pour déterminer ses vaps
  int m = 50; // nombre d'itération de l'algorithme QR 
  schur(C, T, Q, N-2, m);
  aff_T(T, N-2);
  vaps(V, T, N-2);
  

  // On en déduit les vecteurs d'onde k possibles
  double *k = (double*)malloc((N-2)*sizeof(double));
  for (i=0; i<(N-2); i++)
  {
    k[i] = sqrt(-(V[i]));
  }
  tri_croissant(k, N-2); 
  aff_k(k, N-2, l);
  // Pour chaque valeur propre k_i le vecteur propre associé est la i-eme colonne de la matrice Q de la décomposition de Schur

  // exemple d'affichage du vep (on est censé retrouver l'enveloppe sous forme d'un sinus, cf solution analytique du probléme avec conditions de Dirichlet aux bords en 1D)
  //exemple pour le mode 3
  
  for (i=0; i<N-2; i++)
    {
      E[i] = Q[i][MODE-1]; // contient le vep associé au 3e mode 
    }

  fstream vep;
  vep.open("vep.txt", ios::out);
  for (i=0; i<N-2; i++)
    {
      vep << E[i] << endl;
    }
  vep.close();


  fstream fgrille;
  fgrille.open("grille.txt", ios::out); 
  for(i=1; i<(N-1); i++)
  {
   fgrille << grille[i] << endl; 
  }
  fgrille.close(); 


  fstream fdata; 
  fdata.open("donnees.txt", ios::out);
  fdata << l << endl; 
  fdata << N << endl; 
  fdata << MODE << endl; 

  // ********** Affichage des résultats - Comparaison avec les solutions analytiques **********
  
  ostringstream pyth; 
    pyth 
    << "import numpy as np\n"
    << "D = loadtxt('donnees.txt')\n"
    << "E = loadtxt('vep.txt')\n"
    << "G = loadtxt('grille.txt')\n"
    << "m = D[2]\n"
    << "l = D[0]\n"
    << "Z = np.linspace(0, l, 1000)\n"
    << "Y = [np.sin((m*np.pi*z)/l) for z in Z]\n"
    << "for i in range(len(E)): \n"
    << "  plot(G[i], E[i], '+')\n"
    << "plot(Z,Y)\n"
    ; 

  make_plot_py(pyth);



  // ********** Libération de la mémoire allouée et fin du programme **********
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
