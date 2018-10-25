#include "mex.h"
#include <stdio.h>
#include <math.h>
#include <string.h>
#define EPS1 1.0e-12
#define EPS2 1.0e-12
#define EPS3 1.0e-12
#define EPS4 1.0e-12
#define infy 10000

#define MPI 3.1415926536
#define MPIMUL2 6.2831853072
#define MPIDIV2 1.5707963268

float RADCURV;

/***********************************************************/
double mod2pi(angle)
   double angle;
{
   while (angle < 0.0) angle = angle + MPIMUL2;
   while (angle >= MPIMUL2) angle = angle - MPIMUL2;
   return angle;
}
/***********************************************************/
int fct_curve(ty,orientation,val,x1,y1,t1,delta,pathx,pathy,patht,n)
   int ty,orientation;
   double val;
   double *x1,*y1,*t1;
   double delta;
   double *pathx,*pathy,*patht;
   int n;
{
   int i;
   double va1,va2,newval,incrt,remain;
   double center_x,center_y;
   double x2,y2,t2;
   int nnew;

   if (ty == 3)
      if (fabs(val/RADCURV)<EPS4) return(0);
   else
      if (fabs(val)<EPS4) return(0);

   switch(ty)
     {
      case 1 : /* circular arc toward the right */
        center_x = *x1 + RADCURV*sin(*t1);
        center_y = *y1 - RADCURV*cos(*t1);
        va1 = *t1+MPIDIV2;
        if (orientation == 1) va2 = va1-val;
                         else va2 = va1+val;
        x2 = center_x + RADCURV*cos(va2);
        y2 = center_y + RADCURV*sin(va2);
        t2 = *t1 - orientation*val;
        
        nnew = val/delta;
        remain = val - nnew*delta;
        nnew = nnew+n;
	
        if (orientation == -1) delta = -delta;
        incrt = 0;
        for (i = n; i<nnew; i++)
	  {
           va1 = va1-delta;
           *(pathx+i) = center_x + RADCURV*cos(va1);
           *(pathy+i) = center_y + RADCURV*sin(va1);
           incrt = incrt - delta;
           *(patht+i) = mod2pi(*t1 + incrt);
	  }
        n = nnew;
        if (remain > fabs(delta)/5.)
	  {
           *(pathx+nnew) = x2;
           *(pathy+nnew) = y2;
           *(patht+nnew) = mod2pi(t2);
           n++;
	  }
        else
	  {
           *(pathx+nnew-1) = x2;
           *(pathy+nnew-1) = y2;
           *(patht+nnew-1) = mod2pi(t2);
	  }
        break;

      case 2 : /* circular arc toward the left */
        center_x = *x1 - RADCURV*sin(*t1);
        center_y = *y1 + RADCURV*cos(*t1);
        va1 = *t1-MPIDIV2;
        if (orientation == 1) va2 = va1+val;
                         else va2 = va1-val;
        x2 = center_x + RADCURV*cos(va2);
        y2 = center_y + RADCURV*sin(va2);
        t2 = *t1 + orientation*val;
        
        nnew = val/delta;
        remain = val - nnew*delta;
        nnew = nnew+n;

        if (orientation == -1) delta = -delta;
        incrt = 0;
        for (i = n; i<nnew; i++)
	  {
           va1 = va1+delta;
           *(pathx+i) = center_x + RADCURV*cos(va1);
           *(pathy+i) = center_y + RADCURV*sin(va1);
           incrt = incrt + delta;
           *(patht+i) = mod2pi(*t1 + incrt);
	  }
        n = nnew;
        if (remain > fabs(delta)/5.)
	  {
           *(pathx+nnew) = x2;
           *(pathy+nnew) = y2;
           *(patht+nnew) = mod2pi(t2);
           n++;
	 }
        else
	  {
           *(pathx+nnew-1) = x2;
           *(pathy+nnew-1) = y2;
           *(patht+nnew-1) = mod2pi(t2);
	  }
        break;

      case 3 : /* straight line */
        x2 = *x1 + orientation*val*cos(*t1);
        y2 = *y1 + orientation*val*sin(*t1);
        *t1 = mod2pi(*t1);
        t2 = *t1;

        va1 = sqrt((x2-*x1)*(x2-*x1)+(y2-*y1)*(y2-*y1));
        i = va1/0.2;
        remain = va1 - i*0.2;
        nnew = n+i;
        newval = 0.2;
        va1 = orientation*cos(*t1);
        va2 = orientation*sin(*t1);
        for (i = n; i<nnew; i++)
	  {
           *(pathx+i) = *x1 + va1*newval;
           *(pathy+i) = *y1 + va2*newval;
           *(patht+i) = *t1;
           newval = newval + 0.2;
	  }
        if (remain > 0.07)
	  {
            *(pathx+nnew) = x2;
            *(pathy+nnew) = y2;
            *(patht+nnew) = t2;
            n = nnew+1;
	  }
        else
	  {
            *(pathx+nnew-1) = x2;
            *(pathy+nnew-1) = y2;
            *(patht+nnew-1) = t2;
            n = nnew;
	  }
     }

   *x1 = x2;
   *y1 = y2;
   *t1 = t2;

   return(n);
}


/***********************************************************/
int constRS(num,t,u,v,x1,y1,t1,delta,pathx,pathy,patht)
   int num;
   double t,u,v;
   double x1,y1,t1,delta;
   double *pathx,*pathy,*patht;
{
   int left,right,straight,fwd,bwd;
   int n;

   *pathx = x1;
   *pathy = y1;
   *patht = t1;
   n = 1;

   right = 1; left = 2; straight = 3;
   fwd = 1; bwd = -1;
   switch(num)
     {

/*   C | C | C   */

       case 1 : 
         n = fct_curve(left,fwd,t,&x1,&y1,&t1,delta,pathx,pathy,patht,1);
         n = fct_curve(right,bwd,u,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         n = fct_curve(left,fwd,v,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         break;

       case 2 : 
         n = fct_curve(left,bwd,t,&x1,&y1,&t1,delta,pathx,pathy,patht,1);
         n = fct_curve(right,fwd,u,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         n = fct_curve(left,bwd,v,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         break;

       case 3 : 
         n = fct_curve(right,fwd,t,&x1,&y1,&t1,delta,pathx,pathy,patht,1);
         n = fct_curve(left,bwd,u,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         n = fct_curve(right,fwd,v,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         break;

       case 4 : 
         n = fct_curve(right,bwd,t,&x1,&y1,&t1,delta,pathx,pathy,patht,1);
         n = fct_curve(left,fwd,u,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         n = fct_curve(right,bwd,v,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         break;

/*   C | C C   */

       case 5 : 
         n = fct_curve(left,fwd,t,&x1,&y1,&t1,delta,pathx,pathy,patht,1);
         n = fct_curve(right,bwd,u,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         n = fct_curve(left,bwd,v,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         break;

       case 6 : 
         n = fct_curve(left,bwd,t,&x1,&y1,&t1,delta,pathx,pathy,patht,1);
         n = fct_curve(right,fwd,u,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         n = fct_curve(left,fwd,v,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         break;

       case 7 : 
         n = fct_curve(right,fwd,t,&x1,&y1,&t1,delta,pathx,pathy,patht,1);
         n = fct_curve(left,bwd,u,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         n = fct_curve(right,bwd,v,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         break;

       case 8 : 
         n = fct_curve(right,bwd,t,&x1,&y1,&t1,delta,pathx,pathy,patht,1);
         n = fct_curve(left,fwd,u,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         n = fct_curve(right,fwd,v,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         break;

/*   C S C   */

       case 9 : 
         n = fct_curve(left,fwd,t,&x1,&y1,&t1,delta,pathx,pathy,patht,1);
         n = fct_curve(straight,fwd,u,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         n = fct_curve(left,fwd,v,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         break;

       case 10 : 
         n = fct_curve(right,fwd,t,&x1,&y1,&t1,delta,pathx,pathy,patht,1);
         n = fct_curve(straight,fwd,u,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         n = fct_curve(right,fwd,v,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         break;

       case 11 : 
         n = fct_curve(left,bwd,t,&x1,&y1,&t1,delta,pathx,pathy,patht,1);
         n = fct_curve(straight,bwd,u,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         n = fct_curve(left,bwd,v,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         break;

       case 12 : 
         n = fct_curve(right,bwd,t,&x1,&y1,&t1,delta,pathx,pathy,patht,1);
         n = fct_curve(straight,bwd,u,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         n = fct_curve(right,bwd,v,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         break;

       case 13 : 
         n = fct_curve(left,fwd,t,&x1,&y1,&t1,delta,pathx,pathy,patht,1);
         n = fct_curve(straight,fwd,u,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         n = fct_curve(right,fwd,v,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         break;

       case 14 : 
         n = fct_curve(right,fwd,t,&x1,&y1,&t1,delta,pathx,pathy,patht,1);
         n = fct_curve(straight,fwd,u,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         n = fct_curve(left,fwd,v,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         break;

       case 15 : 
         n = fct_curve(left,bwd,t,&x1,&y1,&t1,delta,pathx,pathy,patht,1);
         n = fct_curve(straight,bwd,u,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         n = fct_curve(right,bwd,v,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         break;

       case 16 : 
         n = fct_curve(right,bwd,t,&x1,&y1,&t1,delta,pathx,pathy,patht,1);
         n = fct_curve(straight,bwd,u,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         n = fct_curve(left,bwd,v,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         break;

/*   C Cu | Cu C   */

       case 17 : 
         n = fct_curve(left,fwd,t,&x1,&y1,&t1,delta,pathx,pathy,patht,1);
         n = fct_curve(right,fwd,u,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         n = fct_curve(left,bwd,u,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         n = fct_curve(right,bwd,v,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         break;

       case 18 : 
         n = fct_curve(right,fwd,t,&x1,&y1,&t1,delta,pathx,pathy,patht,1);
         n = fct_curve(left,fwd,u,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         n = fct_curve(right,bwd,u,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         n = fct_curve(left,bwd,v,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         break;

       case 19 : 
         n = fct_curve(left,bwd,t,&x1,&y1,&t1,delta,pathx,pathy,patht,1);
         n = fct_curve(right,bwd,u,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         n = fct_curve(left,fwd,u,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         n = fct_curve(right,fwd,v,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         break;

       case 20 : 
         n = fct_curve(right,bwd,t,&x1,&y1,&t1,delta,pathx,pathy,patht,1);
         n = fct_curve(left,bwd,u,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         n = fct_curve(right,fwd,u,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         n = fct_curve(left,fwd,v,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         break;

/*   C | Cu Cu | C   */

       case 21 : 
         n = fct_curve(left,fwd,t,&x1,&y1,&t1,delta,pathx,pathy,patht,1);
         n = fct_curve(right,bwd,u,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         n = fct_curve(left,bwd,u,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         n = fct_curve(right,fwd,v,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         break;

       case 22 : 
         n = fct_curve(right,fwd,t,&x1,&y1,&t1,delta,pathx,pathy,patht,1);
         n = fct_curve(left,bwd,u,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         n = fct_curve(right,bwd,u,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         n = fct_curve(left,fwd,v,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         break;

       case 23 : 
         n = fct_curve(left,bwd,t,&x1,&y1,&t1,delta,pathx,pathy,patht,1);
         n = fct_curve(right,fwd,u,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         n = fct_curve(left,fwd,u,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         n = fct_curve(right,bwd,v,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         break;

       case 24 : 
         n = fct_curve(right,bwd,t,&x1,&y1,&t1,delta,pathx,pathy,patht,1);
         n = fct_curve(left,fwd,u,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         n = fct_curve(right,fwd,u,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         n = fct_curve(left,bwd,v,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         break;

/*   C | C2 S C   */

       case 25 : 
         n = fct_curve(left,fwd,t,&x1,&y1,&t1,delta,pathx,pathy,patht,1);
         n = fct_curve(right,bwd,MPIDIV2,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         n = fct_curve(straight,bwd,u,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         n = fct_curve(left,bwd,v,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         break;

       case 26 : 
         n = fct_curve(right,fwd,t,&x1,&y1,&t1,delta,pathx,pathy,patht,1);
         n = fct_curve(left,bwd,MPIDIV2,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         n = fct_curve(straight,bwd,u,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         n = fct_curve(right,bwd,v,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         break;

       case 27 : 
         n = fct_curve(left,bwd,t,&x1,&y1,&t1,delta,pathx,pathy,patht,1);
         n = fct_curve(right,fwd,MPIDIV2,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         n = fct_curve(straight,fwd,u,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         n = fct_curve(left,fwd,v,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         break;

       case 28 : 
         n = fct_curve(right,bwd,t,&x1,&y1,&t1,delta,pathx,pathy,patht,1);
         n = fct_curve(left,fwd,MPIDIV2,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         n = fct_curve(straight,fwd,u,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         n = fct_curve(right,fwd,v,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         break;

       case 29 : 
         n = fct_curve(left,fwd,t,&x1,&y1,&t1,delta,pathx,pathy,patht,1);
         n = fct_curve(right,bwd,MPIDIV2,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         n = fct_curve(straight,bwd,u,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         n = fct_curve(right,bwd,v,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         break;

       case 30 : 
         n = fct_curve(right,fwd,t,&x1,&y1,&t1,delta,pathx,pathy,patht,1);
         n = fct_curve(left,bwd,MPIDIV2,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         n = fct_curve(straight,bwd,u,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         n = fct_curve(left,bwd,v,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         break;

       case 31 : 
         n = fct_curve(left,bwd,t,&x1,&y1,&t1,delta,pathx,pathy,patht,1);
         n = fct_curve(right,fwd,MPIDIV2,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         n = fct_curve(straight,fwd,u,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         n = fct_curve(right,fwd,v,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         break;

       case 32 : 
         n = fct_curve(right,bwd,t,&x1,&y1,&t1,delta,pathx,pathy,patht,1);
         n = fct_curve(left,fwd,MPIDIV2,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         n = fct_curve(straight,fwd,u,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         n = fct_curve(left,fwd,v,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         break;

/*   C | C2 S C2 | C   */

       case 33 : 
         n = fct_curve(left,fwd,t,&x1,&y1,&t1,delta,pathx,pathy,patht,1);
         n = fct_curve(right,bwd,MPIDIV2,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         n = fct_curve(straight,bwd,u,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         n = fct_curve(left,bwd,MPIDIV2,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         n = fct_curve(right,fwd,v,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         break;

       case 34 : 
         n = fct_curve(right,fwd,t,&x1,&y1,&t1,delta,pathx,pathy,patht,1);
         n = fct_curve(left,bwd,MPIDIV2,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         n = fct_curve(straight,bwd,u,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         n = fct_curve(right,bwd,MPIDIV2,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         n = fct_curve(left,fwd,v,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         break;

       case 35 : 
         n = fct_curve(left,bwd,t,&x1,&y1,&t1,delta,pathx,pathy,patht,1);
         n = fct_curve(right,fwd,MPIDIV2,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         n = fct_curve(straight,fwd,u,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         n = fct_curve(left,fwd,MPIDIV2,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         n = fct_curve(right,bwd,v,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         break;

       case 36 : 
         n = fct_curve(right,bwd,t,&x1,&y1,&t1,delta,pathx,pathy,patht,1);
         n = fct_curve(left,fwd,MPIDIV2,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         n = fct_curve(straight,fwd,u,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         n = fct_curve(right,fwd,MPIDIV2,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         n = fct_curve(left,bwd,v,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         break;

/*   C C | C   */

       case 37 : 
         n = fct_curve(left,fwd,t,&x1,&y1,&t1,delta,pathx,pathy,patht,1);
         n = fct_curve(right,fwd,u,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         n = fct_curve(left,bwd,v,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         break;

       case 38 : 
         n = fct_curve(right,fwd,t,&x1,&y1,&t1,delta,pathx,pathy,patht,1);
         n = fct_curve(left,fwd,u,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         n = fct_curve(right,bwd,v,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         break;

       case 39 : 
         n = fct_curve(left,bwd,t,&x1,&y1,&t1,delta,pathx,pathy,patht,1);
         n = fct_curve(right,bwd,u,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         n = fct_curve(left,fwd,v,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         break;

       case 40 : 
         n = fct_curve(right,bwd,t,&x1,&y1,&t1,delta,pathx,pathy,patht,1);
         n = fct_curve(left,bwd,u,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         n = fct_curve(right,fwd,v,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         break;

/*   C S C2 | C   */

       case 41 : 
         n = fct_curve(left,fwd,t,&x1,&y1,&t1,delta,pathx,pathy,patht,1);
         n = fct_curve(straight,fwd,u,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         n = fct_curve(right,fwd,MPIDIV2,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         n = fct_curve(left,bwd,v,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         break;

       case 42 : 
         n = fct_curve(right,fwd,t,&x1,&y1,&t1,delta,pathx,pathy,patht,1);
         n = fct_curve(straight,fwd,u,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         n = fct_curve(left,fwd,MPIDIV2,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         n = fct_curve(right,bwd,v,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         break;

       case 43 : 
         n = fct_curve(left,bwd,t,&x1,&y1,&t1,delta,pathx,pathy,patht,1);
         n = fct_curve(straight,bwd,u,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         n = fct_curve(right,bwd,MPIDIV2,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         n = fct_curve(left,fwd,v,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         break;

       case 44 : 
         n = fct_curve(right,bwd,t,&x1,&y1,&t1,delta,pathx,pathy,patht,1);
         n = fct_curve(straight,bwd,u,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         n = fct_curve(left,bwd,MPIDIV2,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         n = fct_curve(right,fwd,v,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         break;

       case 45 : 
         n = fct_curve(left,fwd,t,&x1,&y1,&t1,delta,pathx,pathy,patht,1);
         n = fct_curve(straight,fwd,u,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         n = fct_curve(left,fwd,MPIDIV2,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         n = fct_curve(right,bwd,v,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         break;

       case 46 : 
         n = fct_curve(right,fwd,t,&x1,&y1,&t1,delta,pathx,pathy,patht,1);
         n = fct_curve(straight,fwd,u,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         n = fct_curve(right,fwd,MPIDIV2,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         n = fct_curve(left,bwd,v,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         break;

       case 47 : 
         n = fct_curve(left,bwd,t,&x1,&y1,&t1,delta,pathx,pathy,patht,1);
         n = fct_curve(straight,bwd,u,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         n = fct_curve(left,bwd,MPIDIV2,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         n = fct_curve(right,fwd,v,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         break;

       case 48 : 
         n = fct_curve(right,bwd,t,&x1,&y1,&t1,delta,pathx,pathy,patht,1);
         n = fct_curve(straight,bwd,u,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         n = fct_curve(right,bwd,MPIDIV2,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         n = fct_curve(left,fwd,v,&x1,&y1,&t1,delta,pathx,pathy,patht,n);
         break;



       default:
         printf("Error: RS curve type %d unknown\n",num);
     }

   return n;
}
   
void mexFunction( int nlhs, mxArray *plhs[],
       int nrhs, const mxArray *prhs[] )
{
   if(nrhs!=10) {
       mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","10 inputs required.");}
   if(nlhs!=3) {
       mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","3 output required.");}
   double t = mxGetScalar(prhs[0]);
   double u = mxGetScalar(prhs[1]);
   double v = mxGetScalar(prhs[2]);
   int num = mxGetScalar(prhs[3]);
   double x1 = mxGetScalar(prhs[4]);
   double y1 = mxGetScalar(prhs[5]);
   double t1 = mxGetScalar(prhs[6]);
   RADCURV   = mxGetScalar(prhs[7]);
   double delta = mxGetScalar(prhs[8]);
   int n  = mxGetScalar(prhs[9]);
   double *outx,*outy,*outt;
   
   double pathx[n], pathy[n], patht[n];

   n  = constRS(num,t,u,v,x1,y1,t1,delta,pathx,pathy,patht);

   plhs[0] = mxCreateDoubleMatrix(1,(mwSize)n,mxREAL);
   outx = mxGetPr(plhs[0]);
   plhs[1] = mxCreateDoubleMatrix(1,(mwSize)n,mxREAL);
   outy = mxGetPr(plhs[1]);
   plhs[2] = mxCreateDoubleMatrix(1,(mwSize)n,mxREAL);
   outt = mxGetPr(plhs[2]);
   memcpy(outx,pathx, n*8);
   memcpy(outy,pathy, n*8);
   memcpy(outt,patht, n*8);

}