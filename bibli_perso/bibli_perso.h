#ifndef BIBLI_PERSO_H
#define BIBLI_PERSO_H


#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <time.h>
#include "bibli_fonctions.h"

using namespace std;


void cholesky(double **,int,  double **, double **); 
void produit_scalaire(double *, double *, double *, int); 
void norme(double *, double *, int); 
void vecteur_canonique(double *, int); 
void copie_mat(double **, double **, int); 
void reflexion(double *, double **, int , int); 
void householder(double **, double **, double **, int);
void schur(double **, double **, double**, int , int); 


void affichage(double *, double*, double**, double **); 
void aff_cholesky(double **, double **); 
void aff_c(double **, int); 
void aff_T(double **, int); 
void aff_k(double *, int);


#endif