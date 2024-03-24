#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <time.h>
#include "bibli_fonctions.h" // A decommenter lorsque je veux faire des plots Python sur les ordis du magistère (Python.h pas reconnu par VSC)
//#include"bibli_perso.h"

using namespace std;

const int N = 20; // nombre de points de la grille
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


void aff_c(double **C, int n)
{
  fstream fc;
  fc.open("c.dat", ios::out);

  int i, j;
  for (i=0; i<n; i++)
  {
    for(j=0; j<n; j++)
    {
      fc << C[i][j] << " "; 
    }
    fc << endl;
  }
  fc << endl;
  fc.close();
}




void aff_T(double **T, int n)
{
  fstream ft;
  ft.open("T.dat", ios::out);

  int i, j;
  for (i=0; i<n; i++)
  {
    for(j=0; j<n; j++)
    {
      ft << T[i][j] << " "; 
    }
    ft << endl;
  }
  ft << endl;
  ft.close();

}


void aff_k(double *k, int n)
{
  fstream fk;
  fk.open("vaps.txt", ios::out);
  int i;
  for (i=0; i<n; i++)
    {
      fk << ((i+1)*M_PI)/L << "   " << k[i] << endl;
    }
  fk.close();
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

  fstream f6;
  f6.open("matrice_passage.txt", ios::out);
  
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


    for (l=0; l<n  ; l++)
    {
      for(p=0; p<n  ; p++)
      {
        f6<< Q[l][p] << " "; 
      }
      f6<< endl;
    }
    f6<< endl;


  }

  copie_mat(A, T, n); // T a priori trigsup avec suffisamment de précision pour m assez grand 

  mat_inv(Q, TMP, n);
  copie_mat(TMP,Q, n);  // Q est la matrice de passage qui contient les vecteurs propres !

   
  

  f6.close();


  free_(Q_TMP, n);
  free_(R_TMP, n);
  free_(TMP, n);
  f.close(); 

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

  int m = 30; // nombre d'itération de l'algorithme QR 

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
  tri_croissant(k, N-2);
  aff_k(k, N-2);
  // Pour chaque valeur propre k_i le vecteur propre associé est la i-eme colonne de la matrice Q de la décomposition de Schur

  // D'abord on calcul uniquement le vep associé au mode 3 pour avoir un beau plot, ensuite on fera un truc plus élégant
  
  for (i=0; i<N-2; i++)
    {
      E[i] = Q[i][5];
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
