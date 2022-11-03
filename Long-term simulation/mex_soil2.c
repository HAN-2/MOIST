#include "mex.h"
#include <math.h>
#include<string.h>
// X2.rhob  X.theta X2.thetac theta_sat X.T
void mex_soil2(double *theta,double *theta_sat, double *T, double *kH, double *thao, double *thaol, double *Dv,
size_t m1, size_t m2, size_t m3, size_t n1, size_t n2, size_t n3)
{
    size_t i;
    double const rhob=1300; 
    double const thetac= 0.2;
   
    

    for (int i=0; i<=m1;i=i+1)
    {
        kH[i] = (0.65 - 0.78 * rhob / 1000 + 0.6 * pow((rhob / 1000), 2)) + (1.06 * (rhob / 1000)) * theta[i] -((0.65 - 0.78 * rhob / 1000 + 0.6 * pow((rhob / 1000),2)) - (0.03 + 0.1 * pow((rhob / 1000),2))) *exp(-pow(((1 + 2.6 * pow(thetac,(-0.5))) * theta[i]),4));
        thao[i] = pow((theta_sat[i]-theta[i]),(7.0/3.0))/pow((theta_sat[i]),2);
        thaol[i]=pow((theta[i]),(7.0/3.0))/pow((theta_sat[i]),2);
        Dv[i]= 2.12/101325.0*pow(((T[i]+273.15)/273.15),1.88)*thao[i]*(theta_sat[i]-theta[i]);
   }
}
    
    


 //kH = (0.65 - 0.78 * rhob / 1000 + 0.6 * (rhob / 1000) ^ 2) + (1.06 * (rhob / 1000)) * theta -((0.65 - 0.78 * rhob / 1000 + 0.6 * (rhob / 1000) ^ 2) - (0.03 + 0.1 * (rhob / 1000) .^ 2)) *exp(-((1 + 2.6 * thetac .^ (-0.5)) * theta) .^ 4);
    //thao = ((X2.theta_sat-X.theta(:,i-1)).^(7/3))./((X2.theta_sat).^2);
    //thaol=((X.theta).^(7/3))./((X2.theta_sat).^2);%0.66*(theta./thetasat).^(8/3);
    



    

void mexFunction(int nlhs, mxArray *plhs[],
                int nrhs, mxArray const *prhs[])
{
    size_t m1,n1,m2,n2,m3,n3;
    double *theta, *theta_sat, *T, *kH, *thao, *thaol, *Dv;
    
    m1=mxGetM(prhs[0]);//res column
    n1=mxGetN(prhs[0]);//res row
    m2=mxGetM(prhs[1]);//sat column
    n2=mxGetN(prhs[1]);//sat row
    m3=mxGetM(prhs[2]);//sat column
    n3=mxGetN(prhs[2]);//sat row
    

    
  
    theta=mxGetPr(prhs[0]);
    theta_sat=mxGetPr(prhs[1]);
    T=mxGetPr(prhs[2]);
    


   plhs[0]=mxCreateDoubleMatrix((mwSize)m1,(mwSize)n1,mxREAL);
   plhs[1]=mxCreateDoubleMatrix((mwSize)m1,(mwSize)n1,mxREAL);
   plhs[2]=mxCreateDoubleMatrix((mwSize)m1,(mwSize)n1,mxREAL);
   plhs[3]=mxCreateDoubleMatrix((mwSize)m1,(mwSize)n1,mxREAL);

   kH=mxGetPr(plhs[0]);
   thao=mxGetPr(plhs[1]);
   thaol =mxGetPr(plhs[2]);
   Dv =mxGetPr(plhs[3]);

   mex_soil2(theta, theta_sat, T, kH, thao,thaol,Dv, m1, m2,m3,n1,n2,n3);
    
}