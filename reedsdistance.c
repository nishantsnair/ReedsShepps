#include "mex.h"
#include <stdio.h>
#include <math.h>

#define EPS1 1.0e-12
#define EPS2 1.0e-12
#define EPS3 1.0e-12
#define EPS4 1.0e-12
#define infy 10000

#define MPI 3.1415926536
#define MPIMUL2 6.2831853072
#define MPIDIV2 1.5707963268

float RADCURV,RADCURVMUL2,RADCURVMUL4;
float SQRADCURV,SQRADCURVMUL2;
/***********************************************************/
double mod2pi(angle)
   double angle;
{
   while (angle < 0.0) angle = angle + MPIMUL2;
   while (angle >= MPIMUL2) angle = angle - MPIMUL2;
   return angle;
}
/***********************************************************/
double c_c_c(x,y,phi,rs,rc,t,u,v)
   double x,y,phi,rs,rc,*t,*u,*v;
{
   
   double a,b,u1,theta,alpha,length_rs;

   a = x-rs;
   b = y+rc;
   if ((fabs(a)<EPS3) && (fabs(b)<EPS3)) return(INFINITY);
   u1 = sqrt(a*a+b*b);
   if (u1>RADCURVMUL4) return(INFINITY);
   theta = atan2(b,a);
   alpha = acos(u1/RADCURVMUL4);
   *t = mod2pi(MPIDIV2 + alpha + theta); 
   *u = mod2pi(MPI-2*alpha);
   *v = mod2pi(phi-*t-*u);

   length_rs = RADCURV*(*t+*u+*v);
   return(length_rs);
}


/***********************************************************/
double c_cc(x,y,phi,rs,rc,t,u,v)
   double x,y,phi,rs,rc,*t,*u,*v;
{
   
   double a,b,u1,theta,alpha,length_rs;

   a = x-rs;
   b = y+rc;
   if ((fabs(a)<EPS3) && (fabs(b)<EPS3)) return(INFINITY);
   u1 = sqrt(a*a+b*b);
   if (u1>RADCURVMUL4) return(INFINITY);
   theta = atan2(b,a);
   alpha = acos(u1/RADCURVMUL4);
   *t = mod2pi(MPIDIV2 + alpha + theta); 
   *u = mod2pi(MPI-2*alpha);
   *v = mod2pi(*t+*u-phi);

   length_rs = RADCURV*(*t+*u+*v);
   return(length_rs);
}


/***********************************************************/
double csca(x,y,phi,rs,rc,t,u,v)
   double x,y,phi,rs,rc,*t,*u,*v;
{
   
   double a,b,length_rs;

   a = x-rs;
   b = y+rc;
   *t = mod2pi(atan2(b,a));
   *u = sqrt(a*a+b*b);
   *v = mod2pi(phi-*t);

   length_rs = RADCURV*(*t+*v) + *u;
   return(length_rs);
}


/***********************************************************/
double cscb(x,y,phi,rs,rc,t,u,v)
   double x,y,phi,rs,rc,*t,*u,*v;
{
   
   
   double a,b,u1,theta,alpha,length_rs;

   a = x+rs;
   b = y-rc;
   u1 = sqrt(a*a+b*b);
   if (u1 < RADCURVMUL2) return(INFINITY);
   theta = atan2(b,a);
   *u = sqrt(u1*u1 - SQRADCURVMUL2);
   alpha = atan2(RADCURVMUL2,*u);
   *t = mod2pi(theta+alpha);
   *v = mod2pi(*t-phi);

   length_rs = RADCURV*(*t+*v) + *u;
   return(length_rs);
}


/***********************************************************/
double ccu_cuc(x,y,phi,rs,rc,t,u,v)
   double x,y,phi,rs,rc,*t,*u,*v;
{
   
   double a,b,u1,theta,alpha,length_rs;

   a = x+rs;
   b = y-rc;
   if ((fabs(a)<EPS3) && (fabs(b)<EPS3)) return(INFINITY);
   u1 = sqrt(a*a+b*b);
   if (u1 > RADCURVMUL4) return(INFINITY);
   theta = atan2(b,a);
   if (u1>RADCURVMUL2)
     {
      alpha = acos((u1/2-RADCURV)/RADCURVMUL2);
      *t = mod2pi(MPIDIV2+theta-alpha);
      *u = mod2pi(MPI-alpha);
      *v = mod2pi(phi-*t+2*(*u));
     }
   else
     {
      alpha = acos((u1/2+RADCURV)/(RADCURVMUL2));
      *t = mod2pi(MPIDIV2+theta+alpha);
      *u = mod2pi(alpha);
      *v = mod2pi(phi-*t+2*(*u));
     }

   length_rs = RADCURV*(2*(*u)+*t+*v);
   return(length_rs);
}


/***********************************************************/
double c_cucu_c(x,y,phi,rs,rc,t,u,v)
   double x,y,phi,rs,rc,*t,*u,*v;
{
   
   
   double a,b,u1,theta,alpha,length_rs,va1,va2;

   a = x+rs;
   b = y-rc;
   if ((fabs(a)<EPS3) && (fabs(b)<EPS3)) return(INFINITY);
   u1 = sqrt(a*a+b*b);
   if (u1 > 6*RADCURV) return(INFINITY);
   theta = atan2(b,a);
   va1 = (5*SQRADCURV - u1*u1/4)/SQRADCURVMUL2;
   if ((va1 < 0.0) || (va1 > 1.0)) return(INFINITY);
   *u = acos(va1);
   va2 = sin(*u);
   alpha = asin(RADCURVMUL2*va2/u1);
   *t = mod2pi(MPIDIV2+theta+alpha);
   *v = mod2pi(*t-phi);

   length_rs = RADCURV*(2*(*u)+*t+*v);
   return(length_rs);
}


/***********************************************************/
double c_c2sca(x,y,phi,rs,rc,t,u,v)
   double x,y,phi,rs,rc,*t,*u,*v;
{
   
   
   double a,b,u1,theta,alpha,length_rs;

   a = x-rs;
   b = y+rc;
   u1 = sqrt(a*a+b*b);
   if (u1 < RADCURVMUL2) return(INFINITY);
   theta = atan2(b,a);
   *u = sqrt(u1*u1-SQRADCURVMUL2) - RADCURVMUL2;
   if (*u < 0.0) return(INFINITY);
   alpha = atan2(RADCURVMUL2,(*u+RADCURVMUL2));
   *t = mod2pi(MPIDIV2+theta+alpha);
   *v = mod2pi(*t+MPIDIV2-phi);

   length_rs = RADCURV*(*t+MPIDIV2+*v) + *u;
   return(length_rs);
}


/***********************************************************/
double c_c2scb(x,y,phi,rs,rc,t,u,v)
   double x,y,phi,rs,rc,*t,*u,*v;
{
   
   double a,b,u1,theta,length_rs;

   a = x+rs;
   b = y-rc;
   u1 = sqrt(a*a+b*b);
   if (u1 < RADCURVMUL2) return(INFINITY);
   theta = atan2(b,a);
   *t = mod2pi(MPIDIV2+theta);
   *u = u1-RADCURVMUL2;
   *v = mod2pi(phi-*t-MPIDIV2);

   length_rs = RADCURV*(*t+MPIDIV2+*v) + *u;
   return(length_rs);
}


/***********************************************************/
double c_c2sc2_c(x,y,phi,rs,rc,t,u,v)
   double x,y,phi,rs,rc,*t,*u,*v;
{
   
   
   double a,b,u1,theta,alpha,length_rs;

   a = x+rs;
   b = y-rc;
   u1 = sqrt(a*a+b*b);
   if (u1 < RADCURVMUL4) return(INFINITY);
   theta = atan2(b,a);
   *u = sqrt(u1*u1-SQRADCURVMUL2) - RADCURVMUL4;
   if (*u < 0.0) return(INFINITY);
   alpha = atan2(RADCURVMUL2,(*u+RADCURVMUL4));
   *t = mod2pi(MPIDIV2+theta+alpha);
   *v = mod2pi(*t-phi);

   length_rs = RADCURV*(*t+MPI+*v) + *u;
   return(length_rs);
}


/***********************************************************/
double cc_c(x,y,phi,rs,rc,t,u,v)
   double x,y,phi,rs,rc,*t,*u,*v;
{
   
   
   double a,b,u1,theta,alpha,length_rs,va;

   a = x-rs;
   b = y+rc;
   if ((fabs(a)<EPS3) && (fabs(b)<EPS3)) return(INFINITY);
   u1 = sqrt(a*a+b*b);
   if (u1>RADCURVMUL4) return(INFINITY);
   theta = atan2(b,a);
   *u = acos((8*SQRADCURV - u1*u1)/(8*SQRADCURV));
   va = sin(*u);
   if (fabs(va)<0.001) va = 0.0;
   if ((fabs(va)<0.001) && (fabs(u1)<0.001)) return(INFINITY);
   alpha = asin(RADCURVMUL2*va/u1);
   *t = mod2pi(MPIDIV2 - alpha + theta); 
   *v = mod2pi(*t-*u-phi);

   length_rs = RADCURV*(*t+*u+*v);
   return(length_rs);
}


/***********************************************************/
double csc2_ca(x,y,phi,rs,rc,t,u,v)
   double x,y,phi,rs,rc,*t,*u,*v;
{
   double a,b,u1,theta,alpha,length_rs;

   a = x-rs;
   b = y+rc;
   u1 = sqrt(a*a+b*b);
   if (u1 < RADCURVMUL2) return(INFINITY);
   theta = atan2(b,a);
   *u = sqrt(u1*u1-SQRADCURVMUL2) - RADCURVMUL2;
   if (*u < 0.0) return(INFINITY);
   alpha = atan2((*u+RADCURVMUL2),RADCURVMUL2);
   *t = mod2pi(MPIDIV2+theta-alpha);
   *v = mod2pi(*t-MPIDIV2-phi);

   length_rs = RADCURV*(*t+MPIDIV2+*v) + *u;
   return(length_rs);
}


/***********************************************************/
double csc2_cb(x,y,phi,rs,rc,t,u,v)
   double x,y,phi,rs,rc,*t,*u,*v;
{
   double a,b,u1,theta,length_rs;

   a = x+rs;
   b = y-rc;
   u1 = sqrt(a*a+b*b);
   if (u1 < RADCURVMUL2) return(INFINITY);
   theta = atan2(b,a);
   *t = mod2pi(theta);
   *u = u1 - RADCURVMUL2;
   *v = mod2pi(-*t-MPIDIV2+phi);

   length_rs = RADCURV*(*t+MPIDIV2+*v) + *u;
   return(length_rs);
}

/***********************************************************/
double reed_shepp(x1,y1,t1,x2,y2,t2,numero,tr,ur,vr)
   double x1,y1,t1,x2,y2,t2;
   double *tr,*ur,*vr;
   int *numero;
{
   double x,y,phi;
   double t,u,v,tn,un,vn;
   int num;
   double var,vard,theta,alpha,dx,dy,length;
   double sphi,cphi;
   double ap,am,b1,b2;

/* coordinate change */
   dx = x2 - x1;
   dy = y2 - y1;
   theta = atan2(dy,dx);
   alpha = theta - t1;
   vard = sqrt(dx*dx+dy*dy);
   x = cos(alpha)*vard;
   y = sin(alpha)*vard;
   phi = t2 - t1;

   sphi = sin(phi);
   cphi = cos(phi);

   ap = RADCURV*sphi;
   am = -RADCURV*sphi;
   b1 = RADCURV*(cphi-1);
   b2 = RADCURV*(cphi+1);
   
/*   C | C | C   */

   length = c_c_c(x,y,phi,ap,b1,&tn,&un,&vn);
   num = 1;
   t = tn; u = un; v = vn;

   var = c_c_c(-x,y,-phi,am,b1,&tn,&un,&vn);
   if (var < length)
     {
       length = var;
       num = 2;
       t = tn; u = un; v = vn;
     }

   var = c_c_c(x,-y,-phi,am,b1,&tn,&un,&vn);
   if (var < length)
     {
       length = var;
       num = 3;
       t = tn; u = un; v = vn;
     }

   var = c_c_c(-x,-y,phi,ap,b1,&tn,&un,&vn);
   if (var < length)
     {
       length = var;
       num = 4;
       t = tn; u = un; v = vn;
     }

/*   C | C C   */

   var = c_cc(x,y,phi,ap,b1,&tn,&un,&vn);
   if (var < length)
     {
       length = var;
       num = 5;
       t = tn; u = un; v = vn;
     }

   var = c_cc(-x,y,-phi,am,b1,&tn,&un,&vn);
   if (var < length)
     {
       length = var;
       num = 6;
       t = tn; u = un; v = vn;
     }

   var = c_cc(x,-y,-phi,am,b1,&tn,&un,&vn);
   if (var < length)
     {
       length = var;
       num = 7;
       t = tn; u = un; v = vn;
     }

   var = c_cc(-x,-y,phi,ap,b1,&tn,&un,&vn);
   if (var < length)
     {
       length = var;
       num = 8;
       t = tn; u = un; v = vn;
     }

/*   C S C   */

   var = csca(x,y,phi,ap,b1,&tn,&un,&vn);
   if (var < length)
     {
       length = var;
       num = 9;
       t = tn; u = un; v = vn;
     }

   var = csca(x,-y,-phi,am,b1,&tn,&un,&vn);
   if (var < length)
     {
       length = var;
       num = 10;
       t = tn; u = un; v = vn;
     }

   var = csca(-x,y,-phi,am,b1,&tn,&un,&vn);
   if (var < length)
     {
       length = var;
       num = 11;
       t = tn; u = un; v = vn;
     }

   var = csca(-x,-y,phi,ap,b1,&tn,&un,&vn);
   if (var < length)
     {
       length = var;
       num = 12;
       t = tn; u = un; v = vn;
     }

   var = cscb(x,y,phi,ap,b2,&tn,&un,&vn);
   if (var < length)
     {
       length = var;
       num = 13;
       t = tn; u = un; v = vn;
     }

   var = cscb(x,-y,-phi,am,b2,&tn,&un,&vn);
   if (var < length)
     {
       length = var;
       num = 14;
       t = tn; u = un; v = vn;
     }

   var = cscb(-x,y,-phi,am,b2,&tn,&un,&vn);
   if (var < length)
     {
       length = var;
       num = 15;
       t = tn; u = un; v = vn;
     }

   var = cscb(-x,-y,phi,ap,b2,&tn,&un,&vn);
   if (var < length)
     {
       length = var;
       num = 16;
       t = tn; u = un; v = vn;
     }

/*   C Cu | Cu C   */

   var = ccu_cuc(x,y,phi,ap,b2,&tn,&un,&vn);
   if (var < length)
     {
       length = var;
       num = 17;
       t = tn; u = un; v = vn;
     }

   var = ccu_cuc(x,-y,-phi,am,b2,&tn,&un,&vn);
   if (var < length)
     {
       length = var;
       num = 18;
       t = tn; u = un; v = vn;
     }

   var = ccu_cuc(-x,y,-phi,am,b2,&tn,&un,&vn);
   if (var < length)
     {
       length = var;
       num = 19;
       t = tn; u = un; v = vn;
     }

   var = ccu_cuc(-x,-y,phi,ap,b2,&tn,&un,&vn);
   if (var < length)
     {
       length = var;
       num = 20;
       t = tn; u = un; v = vn;
     }

/*   C | Cu Cu | C   */

   var = c_cucu_c(x,y,phi,ap,b2,&tn,&un,&vn);
   if (var < length)
     {
       length = var;
       num = 21;
       t = tn; u = un; v = vn;
     }

   var = c_cucu_c(x,-y,-phi,am,b2,&tn,&un,&vn);
   if (var < length)
     {
       length = var;
       num = 22;
       t = tn; u = un; v = vn;
     }

   var = c_cucu_c(-x,y,-phi,am,b2,&tn,&un,&vn);
   if (var < length)
     {
       length = var;
       num = 23;
       t = tn; u = un; v = vn;
     }

   var = c_cucu_c(-x,-y,phi,ap,b2,&tn,&un,&vn);
   if (var < length)
     {
       length = var;
       num = 24;
       t = tn; u = un; v = vn;
     }

/*   C | C2 S C   */

   var = c_c2sca(x,y,phi,ap,b1,&tn,&un,&vn);
   if (var < length)
     {
       length = var;
       num = 25;
       t = tn; u = un; v = vn;
     }

   var = c_c2sca(x,-y,-phi,am,b1,&tn,&un,&vn);
   if (var < length)
     {
       length = var;
       num = 26;
       t = tn; u = un; v = vn;
     }

   var = c_c2sca(-x,y,-phi,am,b1,&tn,&un,&vn);
   if (var < length)
     {
       length = var;
       num = 27;
       t = tn; u = un; v = vn;
     }

   var = c_c2sca(-x,-y,phi,ap,b1,&tn,&un,&vn);
   if (var < length)
     {
       length = var;
       num = 28;
       t = tn; u = un; v = vn;
     }

   var = c_c2scb(x,y,phi,ap,b2,&tn,&un,&vn);
   if (var < length)
     {
       length = var;
       num = 29;
       t = tn; u = un; v = vn;
     }

   var = c_c2scb(x,-y,-phi,am,b2,&tn,&un,&vn);
   if (var < length)
     {
       length = var;
       num = 30;
       t = tn; u = un; v = vn;
     }

   var = c_c2scb(-x,y,-phi,am,b2,&tn,&un,&vn);
   if (var < length)
     {
       length = var;
       num = 31;
       t = tn; u = un; v = vn;
     }

   var = c_c2scb(-x,-y,phi,ap,b2,&tn,&un,&vn);
   if (var < length)
     {
       length = var;
       num = 32;
       t = tn; u = un; v = vn;
     }

/*   C | C2 S C2 | C   */

   var = c_c2sc2_c(x,y,phi,ap,b2,&tn,&un,&vn);
   if (var < length)
     {
       length = var;
       num = 33;
       t = tn; u = un; v = vn;
     }

   var = c_c2sc2_c(x,-y,-phi,am,b2,&tn,&un,&vn);
   if (var < length)
     {
       length = var;
       num = 34;
       t = tn; u = un; v = vn;
     }

   var = c_c2sc2_c(-x,y,-phi,am,b2,&tn,&un,&vn);
   if (var < length)
     {
       length = var;
       num = 35;
       t = tn; u = un; v = vn;
     }

   var = c_c2sc2_c(-x,-y,phi,ap,b2,&tn,&un,&vn);
   if (var < length)
     {
       length = var;
       num = 36;
       t = tn; u = un; v = vn;
     }

/*   C C | C   */

   var = cc_c(x,y,phi,ap,b1,&tn,&un,&vn);
   if (var < length)
     {
       length = var;
       num = 37;
       t = tn; u = un; v = vn;
     }

   var = cc_c(x,-y,-phi,am,b1,&tn,&un,&vn);
   if (var < length)
     {
       length = var;
       num = 38;
       t = tn; u = un; v = vn;
     }

   var = cc_c(-x,y,-phi,am,b1,&tn,&un,&vn);
   if (var < length)
     {
       length = var;
       num = 39;
       t = tn; u = un; v = vn;
     }

   var = cc_c(-x,-y,phi,ap,b1,&tn,&un,&vn);
   if (var < length)
     {
       length = var;
       num = 40;
       t = tn; u = un; v = vn;
     }

/*   C S C2 | C   */

   var = csc2_ca(x,y,phi,ap,b1,&tn,&un,&vn);
   if (var < length)
     {
       length = var;
       num = 41;
       t = tn; u = un; v = vn;
     }

   var = csc2_ca(x,-y,-phi,am,b1,&tn,&un,&vn);
   if (var < length)
     {
       length = var;
       num = 42;
       t = tn; u = un; v = vn;
     }

   var = csc2_ca(-x,y,-phi,am,b1,&tn,&un,&vn);
   if (var < length)
     {
       length = var;
       num = 43;
       t = tn; u = un; v = vn;
     }

   var = csc2_ca(-x,-y,phi,ap,b1,&tn,&un,&vn);
   if (var < length)
     {
       length = var;
       num = 44;
       t = tn; u = un; v = vn;
     }

   var = csc2_cb(x,y,phi,ap,b2,&tn,&un,&vn);
   if (var < length)
     {
       length = var;
       num = 45;
       t = tn; u = un; v = vn;
     }

   var = csc2_cb(x,-y,-phi,am,b2,&tn,&un,&vn);
   if (var < length)
     {
       length = var;
       num = 46;
       t = tn; u = un; v = vn;
     }

   var = csc2_cb(-x,y,-phi,am,b2,&tn,&un,&vn);
   if (var < length)
     {
       length = var;
       num = 47;
       t = tn; u = un; v = vn;
     }

   var = csc2_cb(-x,-y,phi,ap,b2,&tn,&un,&vn);
   if (var < length)
     {
       length = var;
       num = 48;
       t = tn; u = un; v = vn;
     }






   *tr = t; *ur = u; *vr = v;
   *numero = num;
   return(length);
}



/***********************************************************/
double min_length_rs(x1,y1,t1,x2,y2,t2,numero,t,u,v)
   int *numero;
   double *t,*u,*v;
   double x1,y1,t1,x2,y2,t2;
{
   double length_rs;

   if ((fabs(x1-x2)<EPS1) && (fabs(y1-y2)<EPS1)
       && (fabs(t1-t2)<EPS1))  length_rs = 0.0;
   else length_rs = reed_shepp(x1,y1,t1,x2,y2,t2,numero,t,u,v);

   return(length_rs);
}







/***********************************************************/



void mexFunction( int nlhs, mxArray *plhs[],
             int nrhs, const mxArray *prhs[] )
{   
    if(nrhs!=7) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","Seven inputs required.");}
    if(nlhs!=5) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","five output required.");}
    double x1 = mxGetScalar(prhs[0]);
    double y1 = mxGetScalar(prhs[1]);
    double t1 = mxGetScalar(prhs[2]);
    double x2 = mxGetScalar(prhs[3]);
    double y2 = mxGetScalar(prhs[4]);
    double t2 = mxGetScalar(prhs[5]);
    RADCURV =mxGetScalar(prhs[6]);
    double *out1,*out2,*out3,*out4,*out5;
    plhs[0] = mxCreateDoubleMatrix((mwSize)1, (mwSize)1, mxREAL);
    out1 = mxGetPr(plhs[0]);
    plhs[1] = mxCreateDoubleMatrix((mwSize)1, (mwSize)1, mxREAL);
    out2 = mxGetPr(plhs[1]);
    plhs[2] = mxCreateDoubleMatrix((mwSize)1, (mwSize)1, mxREAL);
    out3 = mxGetPr(plhs[2]);
    plhs[3] = mxCreateDoubleMatrix((mwSize)1, (mwSize)1, mxREAL);
    out4 = mxGetPr(plhs[3]);
    plhs[4] = mxCreateDoubleMatrix((mwSize)1, (mwSize)1, mxREAL);
    out5 = mxGetPr(plhs[4]);
	RADCURVMUL2 = 2 * RADCURV;
	RADCURVMUL4 = 4 * RADCURV;
	SQRADCURV   = RADCURV * RADCURV;
	SQRADCURVMUL2 = 4 * RADCURV * RADCURV;
	double t, u, v;
	int numero;
    *out1=min_length_rs(x1,y1,t1,x2,y2,t2,&numero,&t,&u,&v);
    *out2=t;
    *out3=u;
    *out4=v;
    *out5=numero;
    
}
