#include <stdio.h>
#include <stdlib.h>
#include <math.h>
void gauss(float signoise,float *xlamnoise);
main()
{
	float xnt,hedeg,xnclim,xlamfdeg,vc,vm,tf,xnp,y,xlamf,yd;
	float t,h,s,yold,ydold,xlamdeg,tgo,xlam,xlamd,xnc;
	float ydd;
	int step,pn;
	FILE *fptr;
	fptr=fopen("DATFIL.TXT","w");
	xnt=0.;
	hedeg=-20.;
	xnclim=999999.;
	pn=0;
	xlamfdeg=-30.;
	vc=4000.;
	vm=3000.;
	tf=10.;
	xnp=3.;
	xlamf=xlamfdeg/57.3;
	y=0.;
	yd=-vm*hedeg/57.3;
	t=0.;
	h=.001;
	s=0.;
L10:
	if(t>(tf-.0001))goto L999;
 	yold=y;
	ydold=yd;
	step=1;
	goto L200;
L66:
	step=2;
 	y=y+h*yd;
 	yd=yd+h*ydd;
	t=t+h;
	goto L200;
L55:	
 	y=.5*(yold+y+h*yd);
 	yd=.5*(ydold+yd+h*ydd);
	s=s+h;
	if(s<.09999)goto L10;
	s=0.;
	xlamdeg=xlam*57.3;
	printf("%6.3f %6.3f %6.3f %6.3f \n",t,y,xnc/32.2,xlamdeg);
	fprintf(fptr,"%6.3f %6.3f %6.3f %6.3f \n",t,y,xnc/32.2,xlamdeg);
	goto L10;
L200:
 	tgo=tf-t+.00001;
 	xlam=y/(vc*tgo);
	xlamd=(y+yd*tgo)/(vc*tgo*tgo);
	if(pn==1)
		xnc=xnp*vc*xlamd;
	else
		xnc=4.*vc*xlamd+xnt+2.*vc*(xlam-xlamf)/tgo;
	if(xnc>xnclim)xnc=xnclim;
	if(xnc<-xnclim)xnc=-xnclim;
	ydd=xnt-xnc;
	if (step<2)
		goto L66;
	else
		goto L55;
L999:
	fclose (fptr);
	return 0;
}
