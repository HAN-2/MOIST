#include "mex.h"
#include <math.h>
#include<string.h>

void mex_soil(double *theta_res,double *theta_sat, double *M, double *alp, double *h, double *N, double *scalingfactor, double *ksat, double *l, double *C, double *S, double *k,
size_t m1, size_t m2, size_t n1, size_t n2,size_t m3, size_t m4, size_t n3, size_t n4,size_t m5, size_t m6, size_t n5, size_t n6,size_t m7, size_t n7,size_t m8, size_t n8,size_t m9, size_t n9)
{
    size_t i;
    double const b=3;  
    for (int i=0; i<=m1;i=i+1)
    {
        C[i]=(theta_res[i]-theta_sat[i])*(-M[i]*pow((1+pow((alp[i]*(-h[i])),N[i])),(-M[i]-1)))*(N[i]*pow(alp[i],N[i])*pow((-h[i]),(N[i]-1)))*scalingfactor[i];
        S[i] = scalingfactor[i]/ pow((1+pow(((alp[i]*(-h[i]))),N[i])),M[i]);
        k[i] = ksat[i]*pow((S[i]/scalingfactor[i]),l[i])*pow((1-pow((1-pow((S[i]/scalingfactor[i]),(1/M[i]))),M[i])),2);
   }
}

//C=1.*(theta_res-theta_sat)*(-M*pow((1+pow((alp*(-h)),N)),(-M-1))).*(N*pow((alp),N)*pow((-h),(N-1)))*scalingfactor;

void mexFunction(int nlhs, mxArray *plhs[],
                int nrhs, mxArray const *prhs[])
{
    size_t m1,n1,m2,n2, m3,n3,m4,n4,m5,n5,m6,n6,m7,n7,m8,n8,m9,n9;
    double *theta_res, *theta_sat, *M, *alp, *h, *N,  *scalingfactor, *ksat, *l, *C, *S, *k;
    
    m1=mxGetM(prhs[0]);//res column
    n1=mxGetN(prhs[0]);//res row
    m2=mxGetM(prhs[1]);//sat column
    n2=mxGetN(prhs[1]);//sat row
    m3=mxGetM(prhs[2]);//C column
    n3=mxGetN(prhs[2]);//C row
    m4=mxGetM(prhs[3]);//M column
    n4=mxGetN(prhs[3]);//M row
    m5=mxGetM(prhs[4]);//alp column
    n5=mxGetN(prhs[4]);//alp row
    m6=mxGetM(prhs[5]);//h column
    n6=mxGetN(prhs[5]);//h row
    m7=mxGetM(prhs[6]);//n column
    n7=mxGetN(prhs[6]);//n row
    m8=mxGetM(prhs[7]);//ksat column
    n8=mxGetN(prhs[7]);//ksat row
    m9=mxGetM(prhs[8]);//l column
    n9=mxGetN(prhs[8]);//l row

    
  
    theta_res=mxGetPr(prhs[0]);
    theta_sat=mxGetPr(prhs[1]);
    M        =mxGetPr(prhs[2]);
    alp      =mxGetPr(prhs[3]);
    h        =mxGetPr(prhs[4]);
    N        =mxGetPr(prhs[5]);
    
    scalingfactor=mxGetPr(prhs[6]);
    ksat     =mxGetPr(prhs[7]);
    l        =mxGetPr(prhs[8]);


   plhs[0]=mxCreateDoubleMatrix((mwSize)m1,(mwSize)n1,mxREAL);
   plhs[1]=mxCreateDoubleMatrix((mwSize)m1,(mwSize)n1,mxREAL);
   plhs[2]=mxCreateDoubleMatrix((mwSize)m1,(mwSize)n1,mxREAL);

   C=mxGetPr(plhs[0]);
   S=mxGetPr(plhs[1]);
   k =mxGetPr(plhs[2]);

   mex_soil(theta_res, theta_sat, M, alp, h, N, scalingfactor,ksat,l, C,S,k, m1, m2, m3,m4,m5,m6,m7,m8,m9,n1,n2,n3,n4,n5,n6,n7,n8,n9);
    
}