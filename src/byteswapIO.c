/*  Copyright Rasmus Munk Larsen, Stanford University, 2003 */

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <ctype.h>
#include "maxsize.h"

#define SINGLE 0
#define DOUBLE 1


/* Macros for byteswapping */
#define SWAP(a,b)   { char _tmp; _tmp=(a); (a)=(b); (b)=_tmp; }
#define SWAP2(a) { char *_p = (char *) &(a); SWAP(_p[0],_p[1]); }
#define SWAP3(a) { char *_p = (char *) &(a); SWAP(_p[0],_p[2]); }
#define SWAP4(a) { char *_p = (char *) &(a); \
                   SWAP(_p[0],_p[3]); SWAP(_p[1],_p[2]); }
#define SWAP8(a) { char *_p = (char *) &(a); \
                   SWAP(_p[0],_p[7]); SWAP(_p[1],_p[6]); \
                   SWAP(_p[2],_p[5]); SWAP(_p[3],_p[4]); }
#define SWAPN(a,n) { char __tmp; char *__ptr; int __i; \
                    __ptr = (char *) &(a); \
                     for(__i=0;__i<((n)/2);__i++) \
                       SWAP(__ptr[__i],__ptr[(n)-1-__i]); }
#define SWAPINT(a)    SWAP4((a))
#define SWAPFLOAT(a)  SWAP4((a))
#define SWAPDOUBLE(a) SWAP8((a))
#define BYTESWAP(a)   { switch(sizeof((a))) { \
                        case 0: case 1: break; \
                        case 2: SWAP2((a)); break; \
                        case 3: SWAP3((a)); break; \
                        case 4: SWAP4((a)); break; \
                        case 8: SWAP8((a)); break; \
                        default: SWAPN(a,sizeof((a))); }; }

#define READINT(a)   fread(&(a),sizeof(int),1,fh); \
                      if (*byteswap) SWAPINT((a)); 
#define READINTS(a,n)  { int __i; fread((a),sizeof(int),n,fh); \
                         if (*byteswap) \
                           for(__i=0;__i<n;__i++) SWAPINT((a)[__i]); } 
#define READFLOAT(a)  fread(&(a),sizeof(float),1,fh); \
                      if (*byteswap) SWAPFLOAT((a)); 
#define READFLOATS(a,n)  { int __i; fread((a),sizeof(float),n,fh); \
                         if (*byteswap) \
                           for(__i=0;__i<n;__i++) SWAPFLOAT((a)[__i]); } 
#define READDOUBLE(a) fread(&(a),sizeof(double),1,fh); \
                      if (*byteswap) SWAPDOUBLE((a)); 
#define READDOUBLES(a,n)  { int __i; fread((a),sizeof(double),n,fh); \
                         if (*byteswap) \
                           for(__i=0;__i<n;__i++) SWAPDOUBLE((a)[__i]); } 


#define WRITEINT(a)   if (*byteswap) SWAPINT((a)); \
                      fwrite(&(a),sizeof(int),1,fh); \
                      if (*byteswap) SWAPINT((a)); 
#define WRITEINTS(a,n)  { int __i; \
                         if (*byteswap) \
                           for(__i=0;__i<n;__i++) SWAPINT((a)[__i]); \
                         fwrite((a),sizeof(int),n,fh); \
                         if (*byteswap) \
                           for(__i=0;__i<n;__i++) SWAPINT((a)[__i]); }
#define WRITEFLOAT(a)   if (*byteswap) SWAPFLOAT((a)); \
                      fwrite(&(a),sizeof(float),1,fh); \
                      if (*byteswap) SWAPFLOAT((a)); 
#define WRITEFLOATS(a,n)  { int __i; \
                         if (*byteswap) \
                           for(__i=0;__i<n;__i++) SWAPFLOAT((a)[__i]); \
                         fwrite((a),sizeof(float),n,fh); \
                         if (*byteswap) \
                           for(__i=0;__i<n;__i++) SWAPFLOAT((a)[__i]);  }
#define WRITEDOUBLE(a)   if (*byteswap) SWAPDOUBLE((a)); \
                      fwrite(&(a),sizeof(double),1,fh); \
                      if (*byteswap) SWAPDOUBLE((a)); 
#define WRITEDOUBLES(a,n)  { int __i; \
                         if (*byteswap) \
                           for(__i=0;__i<n;__i++) SWAPDOUBLE((a)[__i]); \
                         fwrite((a),sizeof(double),n,fh); \
                         if (*byteswap) \
                           for(__i=0;__i<n;__i++) SWAPDOUBLE((a)[__i]); } 



void C_readmesh(FILE *fh, int *precision, int *byteswap, double rthetw[],
		int *ldr, int *N_points);
void C_readfg(FILE *fh, int *precision, int *byteswap, double f1[], double f2[],
	      int *ldf, double g1[], double g2[], int *ldg,
	      int *M_kers, int *N_points, int *N_rad, int *M_nl,
	      int *N_theta, int *M_lm, int inl[], int ilm[]);
void C_rdamdl(FILE *fh, int *byteswap, double x[], double aa[], double data[],
	      int *n, int *nmod, int *ivar, int *icry, int *ia); 
void C_rdmatr(FILE *fh, int *byteswap, int *m, int *n, double a[], int *lda);
void C_wrmatr(FILE *fh, int *byteswap, int *m, int *n, double a[], int *lda);
void C_writebidiag(FILE *fh, int *byteswap, int *n_points, int *m_kers, int *n_iter, 
		   double bidiag[],  int *ldbidiag, double hhvec[], int *ldhhvec, 
		   double U[], int *ldu, double V[], int *ldv);
void C_readbidiag(FILE *fh, int *byteswap, int *n_points, int *m_kers, int *n_iter, 
		  double bidiag[],  int *ldbidiag, double hhvec[], int *ldhhvec, 
		  double U[], int *ldu, double V[], int *ldv);
void C_writeavker(FILE *fh, int *byteswap, int *n_targets, int *n_points, 
		   int *icase, double trdoff[], double targetparms[], 
		   int *ldtargetparms, double rot[], double sigma_rot[], 
		  double depart[], double avker[], int *ldavker);


char *f2c_string(char *s, int len)
{
  /* Copy Fortran character* array to zero terminated string. */
  int i; char *buf;

  if (len<1) {
    printf("f2c_string called with len<=0.\n");
    exit(-1);
  }
  buf = (char *) malloc(len + 1);
  if (buf == NULL) {
    fprintf(stderr, "f2c_string: malloc failed\n");
    exit(-1);
  }
  memcpy(buf, s, len);
  buf[len] = '\0';
  /* Trim trailing non-graphical characters (Fortran padding) */
  i = len - 1;
  while (i >= 0 && !isgraph(buf[i])) i--;
  buf[i + 1] = '\0';
  return buf;
}


void c_readmesh_(int *precision, int *byteswap, double rthetw[],
		 int *ldr, int *N_points, char *filename, int len) 
{
  char *cname;
  FILE *fh;  
  
  cname = f2c_string(filename, len);
  if ( (fh = fopen(cname,"rb")) == NULL) {
    fprintf(stderr,"Failed to open file %s: ",cname);
    perror(NULL);
    exit(-1);
  }
  C_readmesh(fh, precision, byteswap, rthetw, ldr, N_points);
  fclose(fh);
  free(cname);
}


void C_readmesh(FILE *fh, int *precision, int *byteswap, double rthetw[],
		int *ldr, int *N_points) 
{
  int dummy;
  int i,j;
  float *buf;

  READINT(dummy);   /* Record marker */
  READINT(*N_points); 
  if (*N_points > MAXPTS) {
    fprintf(stderr,"Too many meshpoints, MAXPTS = %d.\n ",MAXPTS);
    exit(-1);
  }
  READINT(dummy);  /* Read number of columns in rthetw */

  /* Read array */
  if (*precision==SINGLE) {
    buf = (float *)valloc(3*(*N_points)*sizeof(float));
    READFLOATS(buf,*N_points*3);
    for (j=0; j<3; j++) 
      for (i=0; i<*N_points; i++) 
	rthetw[j*(*ldr) + i] = (double) buf[j * (*N_points) + i];
    free(buf);
  }
  else {
    for (j=0; j<3; j++) {
      READDOUBLES(&rthetw[j*(*ldr)],*N_points);
    }        
  }
}


void c_readfg_(int *precision, int *byteswap, double f1[], double f2[],
	       int *ldf, double g1[], double g2[], int *ldg,
	       int *M_kers, int *N_points, int *N_rad, int *M_nl,
	       int *N_theta, int *M_lm, int inl[], int ilm[],
	       char *filename, int len) 
{
  char *cname;
  FILE *fh;

  cname = f2c_string(filename, len);
  if ( (fh = fopen(cname,"rb")) == NULL) {
    fprintf(stderr,"Failed to open file %s: ",cname);
    perror(NULL);
    exit(-1);
  }
  C_readfg(fh, precision, byteswap,  f1,  f2, ldf,  g1,  g2, ldg, M_kers, 
	   N_points, N_rad, M_nl,N_theta, M_lm, inl, ilm);
  fclose(fh);
  free(cname);
}


void C_readfg(FILE *fh, int *precision, int *byteswap, double f1[], 
	      double f2[],int *ldf, double g1[], double g2[], int *ldg,
	      int *M_kers, int *N_points, int *N_rad, int *M_nl,
	      int *N_theta, int *M_lm, int inl[], int ilm[]) 
{  
  int dummy;
  int i,j;
  float *buf;  

  /**** Begin record ****/
  READINT(dummy);   /* Read record marker */
  READINT(*M_kers);
  if (*M_kers > MAXNLM) {
    fprintf(stderr,"Too many modes, MAXNLM = %d.\n ",MAXNLM);
    exit(-1);
  }
  READINT(*N_points); 
  if (*N_points > MAXPTS) {
    fprintf(stderr,"Too many meshpoints, MAXPTS = %d.\n ",MAXPTS);
    exit(-1);
  }
  READINT(dummy);   /* Read record marker */
  /**** End  record ****/


  /**** Begin record ****/
  READINT(dummy);   /* Read record marker */
  READINT(dummy);
  READINTS(inl,*M_kers);
  READINTS(ilm,*M_kers);
  READINT(dummy);   /* Read record marker */
  /**** End record ****/

  /**** Begin record ****/
  READINT(dummy);   /* Read record marker */
  READINT(*N_rad);
  if (*N_rad > MAXRAD) {
    fprintf(stderr,"Too many radial meshpoints, N_rad = %d > MAXRAD = %d.\n ",
	    *N_rad,MAXRAD);
    exit(-1);
  }
  READINT(*M_nl); 
  if (*M_nl > MAXNL) {
    fprintf(stderr,"Too many multiplets, M_nl = %d > MAXNL = %d.\n ",
	    *M_nl,MAXNL);
    exit(-1);
  }
  READINT(dummy);   /* Read record marker */
  /**** End record ****/

  /**** Begin record ****/
  READINT(dummy);   /* Read record marker */
  if (*precision==SINGLE) {
    buf = (float *)valloc((*N_rad)*(*M_nl)*sizeof(float));
    READFLOATS(buf,(*N_rad)*(*M_nl));
    for (j=0; j<*M_nl; j++) 
      for (i=0; i<*N_rad; i++) 
	f1[j*(*ldf) + i] = (double) buf[j * (*N_rad) + i];
    free(buf);
  }
  else {
    for (j=0; j<*M_nl; j++) 
      READFLOATS(&f1[j*(*ldf)],*N_rad);
  }
  READINT(dummy);   /* Read record marker */
  /**** End record ****/

  /**** Begin record ****/
  READINT(dummy);   /* Read record marker */
  if (*precision==SINGLE) {
    buf = (float *)valloc((*N_rad)*(*M_nl)*sizeof(float));
    READFLOATS(buf,(*N_rad)*(*M_nl));
    for (j=0; j<*M_nl; j++) 
      for (i=0; i<*N_rad; i++) 
	f2[j*(*ldf) + i] = (double) buf[j * (*N_rad) + i];
    free(buf);
  }
  else {
    for (j=0; j<*M_nl; j++) 
      READFLOATS(&f2[j*(*ldf)],*N_rad);
  }
  READINT(dummy);   /* Read record marker */
  /**** End record ****/

  /**** Begin record ****/
  READINT(dummy);   /* Read record marker */
  READINT(*N_theta);
  if (*N_theta > MAXTHETA) {
    fprintf(stderr, "Too many meshpoints in latitude, N_theta = %d > MAXTHETA = %d.\n ", *N_theta,MAXTHETA);
    exit(-1);
  }
  READINT(*M_lm); 
  if (*M_lm > MAXLM) {
    fprintf(stderr,"Too many multiplets, M_lm = %d > MAXLM = %d.\n ",
	    *M_lm,MAXLM);
    exit(-1);
  }
  READINT(dummy);   /* Read record marker */
  /**** End record ****/

  /**** Begin record ****/
  READINT(dummy);   /* Read record marker */
  if (*precision==SINGLE) {
    buf = (float *)valloc((*N_theta)*(*M_lm)*sizeof(float));
    READFLOATS(buf,(*N_theta)*(*M_lm));
    for (j=0; j<*M_lm; j++) 
      for (i=0; i<*N_theta; i++) 
	g1[j*(*ldg) + i] = (double) buf[j * (*N_theta) + i];
    free(buf);
  }
  else {
    for (j=0; j<*M_lm; j++) 
      READFLOATS(&g1[j*(*ldg)],*N_theta);
  }
  READINT(dummy);   /* Read record marker */
  /**** End record ****/

  /**** Begin record ****/
  READINT(dummy);   /* Read record marker */
  if (*precision==SINGLE) {
    buf = (float *)valloc((*N_theta)*(*M_lm)*sizeof(float));
    READFLOATS(buf,(*N_theta)*(*M_lm));
    for (j=0; j<*M_lm; j++) 
      for (i=0; i<*N_theta; i++) 
	g2[j*(*ldg) + i] = (double) buf[j * (*N_theta) + i];
    free(buf);
  }
  else {
    for (j=0; j<*M_lm; j++) 
      READFLOATS(&g2[j*(*ldg)],*N_theta);
  }
  READINT(dummy);   /* Read record marker */
  /**** End record ****/
}


void c_rdamdl_(int *byteswap, double x[], double aa[], double data[],
	       int *n, int *nmod, int *ivar, int *icry, int *ia,
	       char *filename, int len) 
{
  char *cname;
  FILE *fh;
  
  cname = f2c_string(filename, len);
  if ( (fh = fopen(cname,"rb")) == NULL) {
    fprintf(stderr,"Failed to open file %s: ",cname);
    perror(NULL);
    exit(-1);
  }
  C_rdamdl(fh,byteswap, x, aa, data, n, nmod, ivar, icry, ia);
  fclose(fh);
  free(cname);
}

void C_rdamdl(FILE *fh, int *byteswap, double x[], double aa[], double data[],
	      int *nn, int *nmod, int *ivar, int *icry, int *ia) 
{  
  int dummy;
  int i;
  double d8;

  READINT(dummy);   /* Record marker */
  READINT(*nmod); 
  READINT(*nn); 
  READDOUBLES(data,8);
  d8 = data[7]+0.1;  
  if (d8 >= 100) {
    *ivar=8;
  }
  else if (d8 >= 10) {
    *ivar = 6;
  }
  else {
    *ivar = 5;
  }
  for (i=0; i<(*nn); i++) {
    READDOUBLE(x[i]);
    READDOUBLES(&aa[i*(*ia)],*ivar);
  }
  READINT(dummy);   /* Record marker */
  *icry = 0;
}



void c_rdmatr_(int *byteswap, int *m, int *n, double a[], int *lda,
	       char *filename, int len) 
{
  char *cname;
  FILE *fh;

  cname = f2c_string(filename, len);
  if ( (fh = fopen(cname,"rb")) == NULL) {
    fprintf(stderr,"Failed to open file %s: ",cname);
    perror(NULL);
    exit(-1);
  }
  C_rdmatr(fh,byteswap, m, n, a, lda);
  fclose(fh);
  free(cname);
}

void C_rdmatr(FILE *fh, int *byteswap, int *m, int *n, double a[], int *lda ) 
{  
  int dummy;
  int i;

  READINT(dummy);   /* Record marker */
  READINT(*m); 
  if ( *lda < *m ) {
    fprintf(stderr,"C_rdmatr: m > lda. Aborting!\n");
    exit(-1);
  }    
  READINT(*n); 
  for (i=0; i<(*n); i++) 
    READDOUBLES(&a[i*(*lda)],*m);
  READINT(dummy);   /* Record marker */
}


void c_wrmatr_(int *byteswap, int *m, int *n, double a[], int *lda,
	       char *filename, int len) 
{
  FILE *fh;
  char *cname;

  cname = f2c_string(filename, len);
  if ( (fh = fopen(cname,"wb")) == NULL) {
    fprintf(stderr,"Failed to open file %s: ",cname);
    perror(NULL);
    exit(-1);
  }
  C_wrmatr(fh, byteswap, m, n, a, lda);
  fclose(fh);
  free(cname);
}

void C_wrmatr(FILE *fh, int *byteswap, int *m, int *n, double a[], int *lda) 
{  
  int reclen;
  int i;

  reclen = 2*sizeof(int) + (*m)*(*n)*sizeof(double);
  WRITEINT(reclen);   /* Record marker */
  WRITEINT(*m); 
  WRITEINT(*n); 
  for (i=0; i<(*n); i++) 
    WRITEDOUBLES(&a[i*(*lda)],*m);
  WRITEINT(reclen);   /* Record marker */
}



void c_readbidiag_(int *byteswap, int *n_points, int *m_kers, int *n_iter, 
		   double bidiag[],  int *ldbidiag, double hhvec[], int *ldhhvec, 
		   double U[], int *ldu, double V[], int *ldv, char *filename, int len) 
{
  char *cname;
  FILE *fh;

  cname = f2c_string(filename, len);
  if ( (fh = fopen(cname,"rb")) == NULL) {
    fprintf(stderr,"Failed to open file %s: ",cname);
    perror(NULL);
    exit(-1);
  }
  C_readbidiag(fh, byteswap, n_points, m_kers, n_iter, bidiag, ldbidiag, 
	       hhvec, ldhhvec, U, ldu, V, ldv);
  fclose(fh);
  free(cname);
}

void C_readbidiag(FILE *fh, int *byteswap, int *n_points, int *m_kers, int *n_iter, 
		  double bidiag[],  int *ldbidiag, double hhvec[], int *ldhhvec, 
		  double U[], int *ldu, double V[], int *ldv) 
{
  int dummy=0,m,n;

  READINT(dummy);   /* Record marker */
  READINT(*n_points); 
  READINT(*m_kers); 
  READINT(*n_iter); 
  READINT(dummy);   /* Record marker */

  m = *n_iter; n = 2;
  C_rdmatr(fh,byteswap,&m,&n,bidiag,ldbidiag);
  m = *m_kers+2; n = 1;
  C_rdmatr(fh,byteswap,&m,&n,hhvec,ldhhvec);
  m = *n_points; n = *n_iter;
  C_rdmatr(fh,byteswap,&m,&n,U,ldu);
  m = *m_kers-1; n = *n_iter;
  C_rdmatr(fh,byteswap,&m,&n,V,ldv);
}
  
 
void c_writebidiag_(int *byteswap, int *n_points, int *m_kers, int *n_iter, 
		    double bidiag[],  int *ldbidiag, double hhvec[], int *ldhhvec, 
		    double U[], int *ldu, double V[], int *ldv, char *filename, int len) 
{
  char *cname;
  FILE *fh;

  cname = f2c_string(filename, len);
  if ( (fh = fopen(cname,"wb")) == NULL) {
    fprintf(stderr,"Failed to open file %s: ",cname);
    perror(NULL);
    exit(-1);
  }
  C_writebidiag(fh, byteswap, n_points, m_kers, n_iter, bidiag,  ldbidiag, hhvec, 
		ldhhvec, U, ldu, V, ldv);
  fclose(fh);
  free(cname);
}

void C_writebidiag(FILE *fh, int *byteswap, int *n_points, int *m_kers, int *n_iter, 
		   double bidiag[],  int *ldbidiag, double hhvec[], int *ldhhvec, 
		   double U[], int *ldu, double V[], int *ldv) 
{
  int reclen,m,n;

  reclen = 3*sizeof(int);
  WRITEINT(reclen);   /* Record marker */
  WRITEINT(*n_points); 
  WRITEINT(*m_kers); 
  WRITEINT(*n_iter); 
  WRITEINT(reclen);   /* Record marker */

  m = *n_iter; n = 2;
  C_wrmatr(fh,byteswap,&m,&n,bidiag,ldbidiag);
  m = *m_kers+2; n = 1;
  C_wrmatr(fh,byteswap,&m,&n,hhvec,ldhhvec);
  m = *n_points; n = *n_iter;
  C_wrmatr(fh,byteswap,&m,&n,U,ldu);
  m = *m_kers-1; n = *n_iter;
  C_wrmatr(fh,byteswap,&m,&n,V,ldv);
}
  
  

void c_writeavker_(int *byteswap, int *n_targets, int *n_points, int *icase,
		   double trdoff[], double targetparms[], int *ldtargetparms,
		   double rot[], double sigma_rot[], double depart[], 
		   double target[], int *ldtarget, char *filename, int len)
{
  char *cname;
  FILE *fh;

  cname = f2c_string(filename, len);
  if ( (fh = fopen(cname,"wb")) == NULL) {
    fprintf(stderr,"Failed to open file %s: ",cname);
    perror(NULL);
    exit(-1);
  }
  C_writeavker(fh, byteswap, n_targets, n_points, icase, trdoff, targetparms, 
	       ldtargetparms, rot, sigma_rot, depart, target, ldtarget);
  fclose(fh);
  free(cname);
}


void c_appendavker_(int *byteswap, int *n_targets, int *n_points, int *icase,
		   double trdoff[], double targetparms[], int *ldtargetparms,
		   double rot[], double sigma_rot[], double depart[], 
		   double target[], int *ldtarget, char *filename, int len)
{
  char *cname;
  FILE *fh;

  cname = f2c_string(filename, len);
  if ( (fh = fopen(cname,"ab")) == NULL) {
    fprintf(stderr,"Failed to open file %s: ",cname);
    perror(NULL);
    exit(-1);
  }
  C_writeavker(fh, byteswap, n_targets, n_points, icase, trdoff, targetparms, 
	       ldtargetparms, rot, sigma_rot, depart, target, ldtarget);
  fclose(fh);
  free(cname);
}



void C_writeavker(FILE *fh, int *byteswap, int *n_targets, int *n_points, 
		   int *icase, double trdoff[], double targetparms[], 
		   int *ldtargetparms, double rot[], double sigma_rot[], 
		   double depart[], double avker[], int *ldavker)
{
  int i,j, reclen;
  double zero = 0.0;

  reclen = 2*sizeof(int) + (10+(*n_points))*sizeof(double);
  for (i=0; i<*n_targets; i++)
  {
    WRITEINT(reclen);   /* Record marker */
    WRITEINT(*icase); 
    WRITEDOUBLE(trdoff[i]); 
    for (j=0; j<4; j++)
    {
      WRITEDOUBLE(targetparms[i+j*(*ldtargetparms)]); 
    }
    WRITEDOUBLE(rot[i]); 
    WRITEDOUBLE(sigma_rot[i]); 
    WRITEDOUBLE(depart[i]); 
    WRITEDOUBLE(zero); 
    WRITEDOUBLE(zero); 
    WRITEINT(*n_points); 
    WRITEDOUBLES(&avker[i*(*ldavker)], *n_points);     
    WRITEINT(reclen);   /* Record marker */
  }
}
