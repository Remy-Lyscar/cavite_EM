#ifndef DEF_AFFICHAGE
#define DEF_AFFICHAGE


#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <time.h>
#include "bibli_fonctions.h"
#include "bibli_perso.h"


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
  fk.open("vaps.dat", ios::out);
  int i;
  for (i=0; i<n; i++)
    {
      fk << k[i] << endl;
    }
  fk.close();
}


#endif