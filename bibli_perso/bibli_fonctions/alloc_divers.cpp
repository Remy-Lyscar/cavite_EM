#include "bibli_fonctions.h"
/*************************************************************************/
double *D_1(
			 int nligne         /* Nombre d'elements */
			 ){
  int i;
  double *mat;
  mat=(double*)malloc(nligne*sizeof(double));
  if(mat==NULL){
    fprintf(stderr,"D_1 mat: mémoire insuffisante\n");
    exit(-1);
  }
  for(i=0;i<nligne;i++) mat[i]=0.;
  return mat;
}
/*************************************************************************/
void f_D_1(
		     double *mat,      /*  Pointeur sur le vecteur */
		     int ligne         /*  Nombre d'elements */
		     ){
  free(mat);
  return;
}
/*************************************************************************/
int *I_1(
		   int nligne          /* Nombre d'elements */
		   ){
  int i;
  int *mat;
  mat=(int*)malloc(nligne*sizeof(int));
  if(mat==NULL){
    fprintf(stderr,"I_1 mat: mémoire insuffisante\n");
    exit(-1);
  }
  for(i=0;i<nligne;i++) mat[i]=0;
  return mat;
}
/*************************************************************************/
void f_I_1(
		  int *mat,      /*  Pointeur sur le vecteur */
		  int ligne      /*  Nombre d'elements */
		  ){
  free(mat);
  return;
}
/*************************************************************************/
char *C_1(
		     int nligne         /* Nombre d'elements */
		     ){
  char *mat;
  mat=(char*)malloc(nligne*sizeof(char));
  if(mat==NULL){
    fprintf(stderr,"C_1 mat: mémoire insuffisante\n");
    exit(-1);
  }
  return mat;
}
/*************************************************************************/
void f_C_1(
		   char *mat,      /*  Pointeur sur le vecteur */
		   int ligne       /*  Nombre d'elements */
		   ){
  free(mat);
  return;
}
/*************************************************************************/
double **D_2(
			  int nligne,         /* Nombre de lignes */
			  int ncolonne        /* Nombre de colonnes */
			  ){
  double **mat;
  int i,j;
  mat=(double**)malloc(nligne*sizeof(double*));
  if(mat==NULL){
    fprintf(stderr,"D_2 mat: mémoire insuffisante\n");
    exit(-1);
  }
  for(i=0;i<nligne;i++){
    mat[i]=(double*)malloc(ncolonne*sizeof(double));
    if(mat[i]==NULL){
      fprintf(stderr,"D_2 mat[%i]: mémoire insuffisante\n",i);
      exit(-1);
    }
  }
  for(i=0;i<nligne;i++) for(j=0;j<ncolonne;j++) mat[i][j]=0.;
  return mat;
}
/*************************************************************************/
void f_D_2(
		       double **mat,      /*  Pointeur sur la matrice */
		       int ligne,         /*  Nombre de lignes */
		       int ncolonne       /*  Nombre de colonnes */
		       ){
  int i;  
  for(i=0;i<ligne;i++)
    free(mat[i]);
  free(mat);
  return;
}
/*************************************************************************/
int **I_2(
		    int nligne,         /* Nombre de lignes */
		    int ncolonne        /* Nombre de colonnes */
		    ){
  int **mat;
  int i,j;
  mat=(int**)malloc(nligne*sizeof(int*));
  if(mat==NULL){
    fprintf(stderr,"I_2 mat: mémoire insuffisante\n");
    exit(-1);
  }
  for(i=0;i<nligne;i++){
    mat[i]=(int*)malloc(ncolonne*sizeof(int));
    if(mat[i]==NULL){
      fprintf(stderr,"I_2 mat[%i]: mémoire insuffisante\n",i);
      exit(-1);
    }
  }
  for(i=0;i<nligne;i++) for(j=0;j<ncolonne;j++) mat[i][j]=0;
  return mat;
}
/*************************************************************************/
void f_I_2(
		    int **mat,         /*  Pointeur sur la matrice */
		    int ligne,         /*  Nombre de lignes */
		    int ncolonne       /*  Nombre de colonnes */
		    ){
  int i;  
  for(i=0;i<ligne;i++)
    free(mat[i]);
  free(mat);
  return;
}
/*************************************************************************/
char **C_2(
		      int nligne,         /* Nombre de lignes */
		      int ncolonne        /* Nombre de colonnes */
		      ){
  char **mat;
  int i;
  mat=(char**)malloc(nligne*sizeof(char*));
  if(mat==NULL){
    fprintf(stderr,"C_2 mat: mémoire insuffisante\n");
    exit(-1);
  }
  for(i=0;i<nligne;i++){
    mat[i]=(char*)malloc(ncolonne*sizeof(char));
    if(mat[i]==NULL){
      fprintf(stderr,"C_2 mat[%i]: mémoire insuffisante\n",i);
      exit(-1);
    }
  }
  return mat;
}
/*************************************************************************/
void f_C_2(
		     char **mat,        /*  Pointeur sur la matrice */
		     int ligne,         /*  Nombre de lignes */
		     int ncolonne       /*  Nombre de colonnes */
		     ){
  int i;  
  for(i=0;i<ligne;i++)
    free(mat[i]);
  free(mat);
  return;
}
/*************************************************************************/
double ***D_3(
			   int n1,         /* Dimension du premier indice */
			   int n2,         /* Dimension du second indice */
			   int n3          /* Dimension du troisieme indice */
			   )
{
  double ***mat;
  int i,j,k;
  mat=(double***)malloc(n1*sizeof(double**));
  if(mat==NULL)
    {
      fprintf(stderr,"D_3 mat: mémoire insuffisante\n");
      exit(-1);
    }
  for(i=0;i<n1;i++)
    {
      mat[i]=(double**)malloc(n2*sizeof(double*));
      if(mat[i]==NULL)
	{
	  fprintf(stderr,"D_3 mat[%i]: mémoire insuffisante\n",i);
	  exit(-1);
	}
      for (j=0;j<n2;j++)
	{
	  mat[i][j]=(double *)malloc(n3*sizeof(double));
	    if(mat[i][j]==NULL)
	      {
		fprintf(stderr,"D_3 mat[%i][%i]: mémoire insuffisante\n",i,j);
		exit(-1);
	      }
	}
    }
  for(i=0;i<n1;i++) for(j=0;j<n2;j++) for(k=0;k<n3;k++) mat[i][j][k]=0.;
  return mat;
}
/*************************************************************************/
void f_D_3(
                       double ***mat,    /* Pointeur sur le tenseur */
                       int n1,           /* Dimension du premier indice */
                       int n2,           /* Dimension du second indice */
                       int n3            /* Dimension du troisieme indice */
		       )
{
  int i,j;  
  for(i=0;i<n1;i++)
    {
      for(j=0;j<n2;j++) free(mat[i][j]);
      free(mat[i]);
    }
  free(mat);
  return;
}
/*************************************************************************/
int ***I_3(
			   int n1,         /* Dimension du premier indice */
			   int n2,         /* Dimension du second indice */
			   int n3          /* Dimension du troisieme indice */
			   )
{
  int ***mat;
  int i,j,k;
  mat=(int***)malloc(n1*sizeof(int**));
  if(mat==NULL)
    {
      fprintf(stderr,"I_3 mat: mémoire insuffisante\n");
      exit(-1);
    }
  for(i=0;i<n1;i++)
    {
      mat[i]=(int**)malloc(n2*sizeof(int*));
      if(mat[i]==NULL)
	{
	  fprintf(stderr,"I_3 mat[%i]: mémoire insuffisante\n",i);
	  exit(-1);
	}
      for (j=0;j<n2;j++)
	{
	  mat[i][j]=(int *)malloc(n3*sizeof(int));
	    if(mat[i][j]==NULL)
	      {
		fprintf(stderr,"I_3 mat[%i][%i]: mémoire insuffisante\n",i,j);
		exit(-1);
	      }
	}
    }
  for(i=0;i<n1;i++) for(j=0;j<n2;j++) for(k=0;k<n3;k++) mat[i][j][k]=0;
  return mat;
}
/*************************************************************************/
void f_I_3(
                       int ***mat,       /* Pointeur sur le tenseur */
                       int n1,           /* Dimension du premier indice */
                       int n2,           /* Dimension du second indice */
                       int n3            /* Dimension du troisieme indice */
		       )
{
  int i,j;  
  for(i=0;i<n1;i++)
    {
      for(j=0;j<n2;j++) free(mat[i][j]);
      free(mat[i]);
    }
  free(mat);
  return;
}
/*************************************************************************/
char ***C_3(
			   int n1,         /* Dimension du premier indice */
			   int n2,         /* Dimension du second indice */
			   int n3          /* Dimension du troisieme indice */
			   )
{
  char ***mat;
  int i,j;
  mat=(char***)malloc(n1*sizeof(char**));
  if(mat==NULL)
    {
      fprintf(stderr,"C_3 mat: mémoire insuffisante\n");
      exit(-1);
    }
  for(i=0;i<n1;i++)
    {
      mat[i]=(char**)malloc(n2*sizeof(char*));
      if(mat[i]==NULL)
	{
	  fprintf(stderr,"C_3 mat[%i]: mémoire insuffisante\n",i);
	  exit(-1);
	}
      for (j=0;j<n2;j++)
	{
	  mat[i][j]=(char *)malloc(n3*sizeof(char));
	    if(mat[i][j]==NULL)
	      {
		fprintf(stderr,"C_3 mat[%i][%i]: mémoire insuffisante\n",i,j);
		exit(-1);
	      }
	}
    }
  return mat;
}
/*************************************************************************/
void f_C_3(
                       char ***mat,      /* Pointeur sur le tenseur */
                       int n1,           /* Dimension du premier indice */
                       int n2,           /* Dimension du second indice */
                       int n3            /* Dimension du troisieme indice */
		       )
{
  int i,j;  
  for(i=0;i<n1;i++)
    {
      for(j=0;j<n2;j++) free(mat[i][j]);
      free(mat[i]);
    }
  free(mat);
  return;
}
/*************************************************************************/
//*************************************************************************/
void ini_I_1(int *x,int ni,...)
{
  int i;
  va_list ap;
  va_start(ap,ni);
  for(i=0;i<ni;i++) x[i]=va_arg(ap,int);
  va_end(ap);
}
//*************************************************************************/
void ini_I_2(int **x,int ni,int nj,...)
{
  int i,j;
  va_list ap;
  va_start(ap,nj);
  for(i=0;i<ni;i++) for(j=0;j<nj;j++) x[i][j]=va_arg(ap,int);
  va_end(ap);
}
//*************************************************************************/
void ini_I_3(int ***x,int ni,int nj,int nk,...)
{
  int i,j,k;
  va_list ap;
  va_start(ap,nk);
  for(i=0;i<ni;i++) for(j=0;j<nj;j++) for(k=0;k<nk;k++) x[i][j][k]=va_arg(ap,int);
  va_end(ap);
}
//*************************************************************************/
//*************************************************************************/
void ini_D_1(double *x,int ni,...)
{
  int i;
  va_list ap;
  va_start(ap,ni);
  for(i=0;i<ni;i++) x[i]=va_arg(ap,double);
  va_end(ap);
}
//*************************************************************************/
void ini_D_2(double **x,int ni,int nj,...)
{
  int i,j;
  va_list ap;
  va_start(ap,nj);
  for(i=0;i<ni;i++) for(j=0;j<nj;j++) x[i][j]=va_arg(ap,double);
  va_end(ap);
}
//*************************************************************************/
void ini_D_3(double ***x,int ni,int nj,int nk,...)
{
  int i,j,k;
  va_list ap;
  va_start(ap,nk);
  for(i=0;i<ni;i++) for(j=0;j<nj;j++) for(k=0;k<nk;k++) x[i][j][k]=va_arg(ap,double);
  va_end(ap);
}
//*************************************************************************/
//*************************************************************************/
void ini_C_1(char *x,int ni,...)
{
  int i;
  va_list ap;
  va_start(ap,ni);
  for(i=0;i<ni;i++) x[i]=va_arg(ap,char);
  va_end(ap);
}
//*************************************************************************/
void ini_C_2(char **x,int ni,int nj,...)
{
  int i,j;
  va_list ap;
  va_start(ap,nj);
  for(i=0;i<ni;i++) for(j=0;j<nj;j++) x[i][j]=va_arg(ap,char);
  va_end(ap);
}
//*************************************************************************/
void ini_C_3(char ***x,int ni,int nj,int nk,...)
{
  int i,j,k;
  va_list ap;
  va_start(ap,nk);
  for(i=0;i<ni;i++) for(j=0;j<nj;j++) for(k=0;k<nk;k++) x[i][j][k]=va_arg(ap,char);
  va_end(ap);
}
//*************************************************************************/
