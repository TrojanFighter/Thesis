#include <stdio.h>
#include <stdlib.h>
#include <math.h>
main()
{
	int order,step;
	float y,yd,yold,ydold,ydd;
	float k0[3],k1[3],k2[3],k3[3];
	float xin,t,s,h,tnew,w,ytheory;
	FILE *fptr;
	fptr=fopen("DATFIL.TXT","w");
	w=20.;
	xin=1.;
	order=2;
	y=0.;
	yd=0.;
	t=s=0.;
	h=.01;
L5:
	if(t>=1.)
		goto L999;
	s=s+h;
	yold=y;
	ydold=yd;
	step=1;
	goto L200;
L40:
	step=2;
	k0[1]=yd;
	k0[2]=ydd;
 	tnew=t+.5*h;
 	y=yold+.5*h*k0[1];
	yd=ydold+.5*h*k0[2];
	goto L200;
L41:
	step=3;
 	k1[1]=yd;
	k1[2]=ydd;
	tnew=t+.5*h;
	y=yold+.5*h*k1[1];
	yd=ydold+.5*h*k1[2];
	goto L200;
L42:
	step=4;
	k2[1]=yd;
	k2[2]=ydd;
	tnew=t+h;
	y=yold+h*k2[1];
	yd=ydold+h*k2[2];
	goto L200;
L43:
	k3[1]=yd;
	k3[2]=ydd;
	t=tnew;
	y=yold+h*(k0[1]+2.*(k1[1]+k2[1])+k3[1])/6.;
	yd=ydold+h*(k0[2]+2.*(k1[2]+k2[2])+k3[2])/6.;
 	if(s>=.01999)
 		goto L5;
 	s=0.;
 	ytheory=(1.-cos(w*t))/w;
 	printf("%10.5f %10.5f %10.5f \n", t,y,ytheory);
	fprintf(fptr,"%10.5f %10.5f %10.5f \n", t,y,ytheory);
	goto L5;
L200:
 	ydd=w*xin-w*w*y;
	if(step==1)
		goto L40;
	if(step==2)
		goto L41;
	if(step==3)
		goto L42;
	else
		goto L43;
L999:
	fclose (fptr);
	return 0;
}
 
 
