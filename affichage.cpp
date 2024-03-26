#ifndef DEF_AFFICHAGE
#define DEF_AFFICHAGE


#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <time.h>
#include "bibli_fonctions.h"
#include "bibli_perso.h"

using namespace std;


void affichage(double *E, double* grille, double**D, double **M, int n)
{
  fstream f; 
  f.open("verif.txt", ios::out); 

  int k;

  f << "Grille unidimensionnelle, il s'agit ici d'une subdivision régulière de la cavité 1D" << endl;
  f<< endl; 
  for(k=0; k<n; k++)
  {
    f << grille[k] << " "; 
  }
  f << endl;
  f << endl;
  f << endl; 

  f << "Vecteur E, pour vérifier l'initialisation" << endl;
  f << endl; 

  for(k=0; k<n-2; k++)
  {
    f << E[k] << " "; 
  }
  f << endl;
  f << endl;
  f << endl;

  f << "Affichage de la matrice de masse" << endl;
  f << endl; 

  int i,j; 
  for (i=0; i<n-2; i++)
  {
    for(j=0; j<n-2; j++)
    {
      f << M[i][j] << " "; 
    }
    f << endl;
  }
  f << endl;
  f << endl;
  f << endl;

  f << "Affichage de la matrice du Laplacien" << endl;
  f << endl;

   for (i=0; i<n-2; i++)
  {
    for(j=0; j<n-2; j++)
    {
      f << D[i][j] << " "; 
    }
    f << endl;
  }
  f << endl;
  f << endl;
  f << endl;

  f.close(); 
}






void aff_cholesky(double **Lt, double **L, int n)
{
  fstream f;
  f.open("cholesky.txt", ios::out);
  
  int i,j; 
  for (i=0; i<n-2; i++)
  {
    for(j=0; j<n-2; j++)
    {
      f << Lt[i][j] << " "; 
    }
    f << endl;
  }
  f << endl;


  for (i=0; i<n-2; i++)
  {
    for(j=0; j<n-2; j++)
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
  fc.open("c.txt", ios::out);

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
  ft.open("T.txt", ios::out);

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


void aff_k(double *k, int n, double l)
{
  fstream fk;
  fk.open("vaps.txt", ios::out);
  fk << "Valeurs de k analytiques" << "       " << "Valeurs de k calculées avec la méthodes des éléments finis" << endl; 
  fk << endl; 
  int i;
  for (i=0; i<n; i++)
    {
      fk << ((i+1)*M_PI)/l << "                                " << k[i] << endl;
    }
  fk.close();
}


#endif
