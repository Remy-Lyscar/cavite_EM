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
