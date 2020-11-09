/*---------------------------------------------------------------------------*/
/* NRUTILS routine, adapted from Numerical Recipes, 2nd ed. */

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <TSpline.h>
#include <algorithm>
#include "mynrutil.hh"
#define NR_END	    1
#define FREE_ARG    char*
#define TOL 1.0e-5

void nrerror(char error_text[]){
    fprintf(stderr, "Numerical Recipes run-time error. . .\n");
    fprintf(stderr, "%s\n", error_text);
    fprintf(stderr, ". . .now exiting to system. . .\n");
    exit(1);
}

/* allocate a float vector with subscript range v[nl..nh] */
float *vector(long nl, long nh){
    float *v;
    
    v = (float *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(float)));
    if (!v) nrerror("allocation failure in vector()");
    return v-nl+NR_END;
}

/* allocate an int vector with subscript range v[nl..nh] */
int *ivector(long nl, long nh){
    int *v;
    
    v = (int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
    if (!v) nrerror("allocation failure in ivector()");
    return v-nl+NR_END;
}

/* allocate a unsigned char vector with subscript range v[nl..nh] */
unsigned char *cvector(long nl, long nh){
  unsigned char *v;
    
  v = (unsigned char *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(unsigned char)));
  if (!v) nrerror("allocation failure in cvector()");
  return v-nl+NR_END;
}

/* allocate a unsigned long vector with subscript range v[nl..nh] */
unsigned long *lvector(long nl, long nh){
    unsigned long *v;
    
    v = (unsigned long *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(long)));
    if (!v) nrerror("allocation failure in lvector()");
    return v-nl+NR_END;
}

/* allocate a double vector with subscript range v[nl..nh] */
double *dvector(long nl, long nh){
    double *v;
    
    v = (double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
    if (!v) nrerror("allocation failure in dvector()");
    return v-nl+NR_END;
}

/* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
float **matrix(long nrl, long nrh, long ncl, long nch){
    long i, nrow=nrh-nrl+1, ncol=nch-ncl+1;
    float **m;
    
    /* allocate pointers to rows */
    m=(float **) malloc((size_t)((nrow+NR_END)*sizeof(float*)));
    if (!m) nrerror("allocation failure 1 in matrix()");
    m += NR_END;
    m -= nrl;
    
    /* allocate rows and set pointers to them */
    m[nrl]=(float *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(float)));
    if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
    m[nrl] += NR_END;
    m[nrl] -= ncl;
    
    for (i=nrl+1;i<=nrh;i++) m[i] = m[i-1] + ncol;
    
    /* return pointer to array of pointers to rows */
    return m;
}

/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
double **dmatrix(long nrl, long nrh, long ncl, long nch){
    long i, nrow=nrh-nrl+1, ncol=nch-ncl+1;
    double **m;
    
    /* allocate pointers to rows */
    m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
    if (!m) nrerror("allocation failure 1 in matrix()");
    m += NR_END;
    m -= nrl;
    
    /* allocate rows and set pointers to them */
    m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
    if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
    m[nrl] += NR_END;
    m[nrl] -= ncl;
    
    for (i=nrl+1;i<=nrh;i++) m[i] = m[i-1] + ncol;
    
    /* return pointer to array of pointers to rows */
    return m;
}

/* allocate an int matrix with subscript range m[nrl..nrh][ncl..nch] */
int **imatrix(long nrl, long nrh, long ncl, long nch){
    long i, nrow=nrh-nrl+1, ncol=nch-ncl+1;
    int **m;
    
    /* allocate pointers to rows */
    m=(int **) malloc((size_t)((nrow+NR_END)*sizeof(int*)));
    if (!m) nrerror("allocation failure 1 in matrix()");
    m += NR_END;
    m -= nrl;
    
    /* allocate rows and set pointers to them */
    m[nrl]=(int *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(int)));
    if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
    m[nrl] += NR_END;
    m[nrl] -= ncl;
    
    for (i=nrl+1;i<=nrh;i++) m[i] = m[i-1] + ncol;
    
    /* return pointer to array of pointers to rows */
    return m;
}

/* point a submatrix [newrl..][newcl..] to a[oldrl..oldrh][oldcl..oldch] */
float **submatrix(float **a, long oldrl, long oldrh, long oldcl, long oldch, 
                  long newrl, long newcl){
    long i, j, nrow=oldrh-oldrl+1, ncol=oldcl-newcl;
    float **m;
    
    /* allocate array of pointers to rows */
    m=(float **) malloc((size_t)((nrow+NR_END)*sizeof(float*)));
    if (!m) nrerror("allocation failure 1 in submatrix()");
    m += NR_END;
    m -= newrl;
    
    /* set pointers to rows */
    for (i=oldrl, j=newrl;i<=oldrh;i++, j++) m[j] = a[i] + ncol;
    
    /* return pointer to array of pointers to rows */
    return m;
}

/* allocate a float matrix m[nrl..nrh][ncl..nch] that points to the matrix
 * declared in the standard C manner as a[nrow][ncol],  where nrow=nrh-nrc+1
 * and ncol=nch-ncl+1.  The routine should be called with the address
 * &a[0][0] as the first argument.
 */
float **convert_matrix(float *a, long nrl, long nrh, long ncl, long nch){
    long i, j, nrow=nrh-nrl+1, ncol=nch-ncl+1;
    float **m;
    
    /* allocate pointers to rows */
    m=(float **) malloc((size_t)((nrow+NR_END)*sizeof(float*)));
    if (!m) nrerror("allocation failure 1 in matrix()");
    m += NR_END;
    m -= nrl;
    
    /* allocate rows and set pointers to them */
    m =(float **) malloc((size_t)((nrow+NR_END)*sizeof(float)));
    if (!m) nrerror("allocation failure 2 in convert_matrix()");
    m += NR_END;
    m -= ncl;
    
    /* set pointers to rows */
    m[nrl] = a - ncl;
    for (i=1, j=nrl+1;i<nrow;i++, j++) m[j] = m[j-1] + ncol;
    
    /* return pointer to array of pointers to rows */
    return m;
}

/* alocate a float 3tensor with range[nrl..nrh][ncl..nch][ndl..ndh] */
float ***f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh){
   long i, j, nrow=nrh-nrl+1, ncol=nch-ncl+1, ndep=ndh-ndl+1;
   float ***t;
    
   /* allocate pointers to pointers to rows */
   t=(float ***) malloc((size_t)((nrow+NR_END)*sizeof(float**)));
   if (!t) nrerror("allocation failure 1 in f3tensor()");
   t += NR_END;
   t -= nrl;
    
   /* allocate pointers to rows and set pointers to them */
   t[nrl]=(float **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(float*)));
   if (!t[nrl]) nrerror("allocation failure 2 in f3tensor()");
   t[nrl] += NR_END;
   t[nrl] -= ncl;
    
   /* allocate rows and set pointers to them */
   t[nrl][ncl]=(float *) malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(float)));
   if (!t[nrl][ncl]) nrerror("allocation failure 3 in f3tensor()");
   t[nrl][ncl] += NR_END;
   t[nrl][ncl] -= ndl;
   
   for(j=ncl+1;j<=nch;j++) t[nrl][j] = t[nrl][j-1]+ndep;
   for(i=nrl+1;i<=nrh;i++) {
       t[i] = t[i-1] + ncol;
       t[i][ncl] = t[i-1][ncl] + ncol*ndep;
       for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
   }
   
   /* return pointer to array of pointers to rows */
   return t;
}

/* free a float vector allocated with vector() */
void free_vector(float *v, long nl, long nh){
    free((FREE_ARG) (v+nl-NR_END));
}

/* free an int vector allocated with ivector() */
void free_ivector(int *v, long nl, long nh){
    free((FREE_ARG) (v+nl-NR_END));
}

/* free an unsigned char vector allocated with cvector() */
void free_cvector(unsigned char *v, long nl, long nh){
    free((FREE_ARG) (v+nl-NR_END));
}

/* free an unsigned long vector allocated with lvector() */
void free_lvector(unsigned long *v, long nl, long nh){
    free((FREE_ARG) (v+nl-NR_END));
}

/* free a double vector allocated with dvector() */
void free_dvector(double *v, long nl, long nh){
    free((FREE_ARG) (v+nl-NR_END));
}

/* free a float matrix allocated by matrix() */
void free_matrix(float **m, long nrl, long nrh, long ncl, long nch){
    free((FREE_ARG) (m[nrl]+ncl-NR_END));
    free((FREE_ARG) (m+nrl-NR_END));
}

/* free a double matrix allocated by dmatrix() */
void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch){
    free((FREE_ARG) (m[nrl]+ncl-NR_END));
    free((FREE_ARG) (m+nrl-NR_END));
}

/* free an int matrix allocated by imatrix() */
void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch){
    free((FREE_ARG) (m[nrl]+ncl-NR_END));
    free((FREE_ARG) (m+nrl-NR_END));
}

/* free a submatrix allocated by submatrix() */
void free_submatrix(float **b, long nrl, long nrh, long ncl, long nch) {
    free((FREE_ARG) (b+nrl-NR_END));
}

/* free a matrix allocated by convert_matrix() */
void free_convert_matrix(float **b, long nrl, long nrh, long ncl, long nch) {
    free((FREE_ARG) (b+nrl-NR_END));
}

/* free a float f3tensor allocated by f3tensor() */
void free_f3tensor(float ***t, long nrl, long nrh, long ncl, long nch, 
                   long ndl, long ndh){
    free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
    free((FREE_ARG) (t[nrl]+ncl-NR_END));
    free((FREE_ARG) (t+nrl-NR_END));
}

/*********************************************************************/
/*  Search a[0:ma-1] by minizing chi^2 for data x[0:ndata-1] and     */
/* y[0:ndata-1] with error sig[0:ndata-1].                           */
/*  The trial function is y=Sum_ia[i]*afunc[i](x).                   */
/* u[0:ndata-1][0:ma-1], v[0:ma-1][0:ma-1], and w[0:ma-1] are used   */
/* for the working space.                                            */
/*  After the excution of this routine, v and w can be used to       */
/* evaluate covariance.                                              */
/*  Parameters a[0:ma-1] and chi-square chisq are calculated in this */
/* routine.                                                          */
/*  User function funcs(x,afunc,ma) evaluates base functions (ma) at */
/* x and store the results into the vector afunc[0:ma-1].            */
/*********************************************************************/

void svdfit(double y[], double sig[], int ndata,
	    double a[], int ma, double **u, double **v, double w[],
	    double *chisq, KTheodata &theo,
	    void(*funcs)(int, double [], int, KTheodata &theo)){

  int j,i;
  double wmax,tmp,thresh,sum,*b,*afunc;

  b=dvector(0,ndata-1);
  afunc=dvector(0,ma-1);

  for(i=0;i<ndata;i++){
    (*funcs)(i,afunc,ma,theo);
    tmp=1.0/sig[i];    
    for(j=0;j<ma;j++)u[i][j]=afunc[j]*tmp;
    b[i]=y[i]*tmp;
  }
  svdcmp(u,ndata,ma,w,v);
  wmax=0.0;
  for(j=0;j<ma;j++)
    if(w[j]>wmax) wmax=w[j];
  thresh=TOL*wmax;
  for(j=0;j<ma;j++)
    if(w[j]<thresh) w[j]=0.0;
  svbksb(u,w,v,ndata,ma,b,a);
  *chisq=0.0;
  for(i=0;i<ndata;i++){
    (*funcs)(i,afunc,ma,theo);
    for(sum=0.0,j=0;j<ma;j++) sum+=a[j]*afunc[j];
    *chisq+=(tmp=(y[i]-sum)/sig[i],tmp*tmp);
  }
  free_dvector(afunc,0,ma-1);
  free_dvector(b,0,ndata-1);
}


/*********************************************************************/
/* Calculate Covariant Matrix cvm[0:ma-1][0:ma-1].                   */
/* Use v[0:ma-1][0:ma-1] and w[0:ma-1] obtained from svdfit.         */
/*********************************************************************/

void svdvar(double **v,int ma,double w[],double **cvm){
  int k,j,i;
  double sum,*wti;

  wti=dvector(0,ma-1);
  for(i=0;i<ma;i++){
    wti[i]=0;
    if(w[i]) wti[i]=1.0/(w[i]*w[i]);
  }
  for(i=0;i<ma;i++){
    for(j=0;j<=i;j++){
      for(sum=0.0,k=0;k<ma;k++) sum+=v[i][k]*v[j][k]*wti[k];
      cvm[j][i]=cvm[i][j]=sum;
    }
  }
  free_dvector(wti,0,ma-1);
}

/**********************************************************************/
/* Solve Ax=b to obtain vector x. A must be SVD decomposed into       */
/* u[0:m-1][0:n-1], w[0:n-1], and v[0:n-1][0:n-1] by calling svdcmp.  */
/* The answer of Ax=b is stored in x[0:n-1].                          */
/**********************************************************************/

void svbksb(double **u, double w[], double **v, int m, int n,
	    double b[], double x[]){
  int jj,j,i;
  double s,*tmp;

  tmp=dvector(0,n-1);
  for(j=0;j<n;j++){ /* Calculate U^tb */
    s=0.0;
    if(w[j]){ /* Only if w[j]!=0, non-zero result is obtained. */
      for(i=0;i<m;i++) s+=u[i][j]*b[i];
      s/=w[j];
    }
    tmp[j]=s;
  }
  for(j=0;j<n;j++){ /* Multiplied by V */
    s=0.0;
    for(jj=0;jj<n;jj++) s+=v[j][jj]*tmp[jj];
    x[j]=s;
  }
  free_dvector(tmp,0,n-1);
}

/* From http://www.cv.cs.ritsumei.ac.jp/~noriko/free/svd.html */

/**************************************************************/
/* Decompose a matrix a[0:m-1][0:n-1] into a=u*w*v^t.         */
/* Memory space **a is used to store the matrix u.            */
/* The matrix v (not v^t) is stored in **v.                   */
/* It must be m>=n. If m<n, a should be a[0:n-1][0:n-1]       */
/* instead of a[0:m-1][0:n-1].                                */
/* The partial matrix a[m:n-1][0:n-1] must be zero matrix.    */
/**************************************************************/

void svdcmp(double **a,int m, int n, double w[], double **v)
{
  double pythag(double a,double b);
  int flag,i,its,j,jj,k,l,nm;
  double anorm,c,f,g,h,s,scale,x,y,z,*rv1;

  /*
  int ii,jjj,aa,bb,cc,ll,eigval[5];
  double vv[5][5],ac[5][5],ww[5],test;*/

  rv1=dvector(0,n-1);
  g=scale=anorm=0.0;
  for (i=0;i<n;i++){
    l=i+1;
    rv1[i]=scale*g;
    g=s=scale=0.0;
    if(i<m){
      for(k=i;k<m;k++) scale +=fabs(a[k][i]);
      if(scale){
        for(k=i;k<m;k++){
          a[k][i]/=scale;
          s+=a[k][i]*a[k][i];
        }
        f=a[i][i];
        g=-SIGN(sqrt(s),f);
        h=f*g-s;
        a[i][i]=f-g;
        for(j=l;j<n;j++){
          for(s=0.0,k=i;k<m;k++) s+= a[k][i]*a[k][j];
          f=s/h;
          for(k=i;k<m;k++) a[k][j]+=f*a[k][i];
        }
        for(k=i;k<m;k++) a[k][i] *=scale;
      }
    }
    w[i]=scale*g;
    g=s=scale=0.0;
    if(i<m && (i != n-1)){
      for(k=l;k<n;k++) scale+=fabs(a[i][k]);
      if(scale){
        for(k=l;k<n;k++){
          a[i][k]/= scale;
          s += a[i][k]*a[i][k];
        }
        f=a[i][l];
        g=-SIGN(sqrt(s),f);
        h=f*g-s;
        a[i][l]=f-g;
        for(k=l;k<n;k++) rv1[k]=a[i][k]/h;
        for(j=l;j<m;j++){
          for (s=0.0,k=l;k<n;k++) s+=a[j][k]*a[i][k];
          for (k=l;k<n;k++) a[j][k] += s*rv1[k];
        }
        for (k=l; k<n; k++) a[i][k] *= scale;
      }
    }
    anorm = FMAX(anorm,(fabs(w[i])+fabs(rv1[i])));
  }
  for(i=n-1;i>=0;i--){
    /******/
    if(i<n-1){
      if(g){
        for(j=l;j<n;j++)
          v[j][i]=(a[i][j]/a[i][l])/g;
        for (j=l;j<n;j++){
          for (s=0.0,k=l;k<n;k++) s += a[i][k]*v[k][j];
          for (k=l;k<n;k++) v[k][j] += s*v[k][i];
        }
      }
      for(j=l;j<n;j++) v[i][j]=v[j][i]=0.0;
    }
    v[i][i]=1.0;
    g=rv1[i];
    l=i;
  }

  for(i=IMIN(m,n)-1;i>=0;i--){
    l=i+1;
    g=w[i];
    for (j=l;j<n;j++) a[i][j]=0.0;
    if (g) {
      g=1.0/g;
      for (j=l; j<n; j++){
        for (s=0.0,k=l;k<m;k++) s += a[k][i]*a[k][j];
        f=(s/a[i][i])*g;
        for (k=i;k<m;k++) a[k][j] += f*a[k][i];
      }
      for(j=i;j<m;j++) a[j][i]*=g;
    }else for (j=i;j<m;j++) a[j][i]=0.0;
    ++a[i][i];
  }
  for (k=n-1;k>=0;k--){
    for (its=1;its<=30;its++){
      flag=1;
      for (l=k;l>=0;l--){
        nm=l-1;
        if ( (double)(fabs(rv1[l])+anorm) == anorm){
          flag=0;
          break;
        }
        if( ( double)(fabs(w[nm])+anorm) == anorm) break;
      }
      if(flag){
        c=0.0;
        s=1.0;
        for (i=l; i<k; i++){
          f=s*rv1[i];
          rv1[i]=c*rv1[i];
          if ((double)(fabs(f)+anorm)== anorm) break;
          g=w[i];
          h=pythag(f,g);
          w[i]=h;
          h=1.0/h;
          c=g*h;
          s=-f*h;
          for(j=0;j<m;j++){
            y=a[j][nm];
            z=a[j][i];
            a[j][nm]=y*c+z*s;
            a[j][i]=z*c-y*s;
          }
        }
      }
      z=w[k];
      if (l==k){
        if(z<0.0){
          w[k]=-z;
          for(j=0;j<n;j++) v[j][k] = -v[j][k];
        }
        break;
      }
      if (its == 30) nrerror("no convergence in 30 svdcmp iterations");
      x=w[l];
      nm=k-1;
      y=w[nm];
      g=rv1[nm];
      h=rv1[k];
      f=( (y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
      g=pythag(f,1.0);
      f=( (x-z)*(x+z)+h*( (y/(f+SIGN(g,f)))-h))/x;
      c=s=1.0;
      for (j=l; j<=nm; j++){
        i=j+1;
        g=rv1[i];
        y=w[i];
        h=s*g;
        g=c*g;
        z=pythag(f,h);
        rv1[j]=z;
        c=f/z;
        s=h/z;
        f=x*c+g*s;
        g=g*c-x*s;
        h=y*s;
        y*=c;
        for (jj=0;jj<n;jj++){
          x=v[jj][j];
          z=v[jj][i];
          v[jj][j]=x*c+z*s;
          v[jj][i]=z*c-x*s;
        }
        z=pythag(f,h);
        w[j]=z;
        if (z) {
          z=1.0/z;
          c=f*z;
          s=h*z;
        }
        f=c*g+s*y;
        x=c*y-s*g;
        for (jj=0;jj<m;jj++){
          y=a[jj][j];
          z=a[jj][i];
          a[jj][j]=y*c+z*s;
          a[jj][i]=z*c-y*s;
        }
      }
      rv1[l]=0.0;
      rv1[k]=f;
      w[k]=x;
    }
  }

  free_dvector(rv1,0,n-1);


}

double pythag(double a,double b)
{
  double absa,absb;
  absa=fabs(a);
  absb=fabs(b);
  if (absa > absb) return absa*sqrt(1.0+SQR(absb/absa));
  else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb)));
}

