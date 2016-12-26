#include <stdio.h>
#include <stdlib.h>
#include <math.h>
main()
{
	float vc,xnt,y,vm,hedeg,tf,xnp,yd,t,h,s,yold,ydold,ydd;
	float tgo,xlamd,xnc;
	int step;
	FILE *fptr;
	fptr=fopen("DATFIL.TXT","w");
	vc=4000.;
	xnt=0.;
	y=0.;
	vm=3000.;
	hedeg=-20.;
	tf=10.;
	xnp=4.;
	yd=-vm*hedeg/57.3;
	t=0.;
	h=.01;
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
	printf("%10.4f %10.4f %10.4f %10.4f\n",t,y,yd,xnc);
	fprintf(fptr,"%10.4f %10.4f %10.4f %10.4f\n",t,y,yd,xnc);
	goto L10;
L200:
 	tgo=tf-t+.00001;
	xlamd=(y+yd*tgo)/(vc*tgo*tgo);
	xnc=xnp*vc*xlamd;
	ydd=xnt-xnc;
	if(step<=1)
		goto L66;
	else
		goto L55;
L999:
 	printf("%10.4f %10.4f %10.4f %10.4f\n",t,y,yd,xnc);
	fclose (fptr);
	return 0;
}
