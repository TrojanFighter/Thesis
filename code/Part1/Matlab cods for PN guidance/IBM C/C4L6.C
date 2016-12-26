#include <stdio.h>
#include <stdlib.h>
#include <math.h>
main()
{
	int step;
	float xnt,xnp,tau,tf,t,s,tp,x1,x2,x3,x4,x5,h;
	float x1old,x2old,x3old,x4old,x5old,x1d,x2d,x3d,x4d,x5d;
	float y1,tgo,xmudnt;
	FILE *fptr;
	fptr=fopen("DATFIL.TXT","w");
	xnt=96.6;
	xnp=3.;
	tau=1.;
	tf=10.;
	t=0.;
	s=0.;
	tp=t+.00001;
	x1=0;
	x2=0;
	x3=1;
	x4=0;
	x5=0.;
	h=.01;
 L10:	if(tp>(tf-.00001))goto L999;
		s=s+h;
		x1old=x1;
		x2old=x2;
		x3old=x3;
		x4old=x4;
		x5old=x5;
		step=1;
		goto L200;
L66:	step=2;
		x1=x1+h*x1d;
		x2=x2+h*x2d;
		x3=x3+h*x3d;
		x4=x4+h*x4d;
		x5=x5+h*x5d;
		tp=tp+h;
		goto L200;
L55:
		x1=(x1old+x1)/2+.5*h*x1d;
		x2=(x2old+x2)/2+.5*h*x2d;
		x3=(x3old+x3)/2+.5*h*x3d;
		x4=(x4old+x4)/2+.5*h*x4d;
		x5=(x5old+x5)/2+.5*h*x5d;
		s=s+h;
		if(s<.099999)goto L10;
		s=0.;
		xmudnt=xnt*sqrt(x5/tgo);
 		printf("%6.3f %6.3f\n",tp,xmudnt);
		fprintf(fptr,"%6.3f %6.3f\n",tp,xmudnt);
		goto L10;
L200:
		x1d=x2;
		x2d=x3;
		y1=(x4-x2)/tau;
		tgo=tp+.00001;
		x3d=xnp*y1/tgo;
		x4d=-y1;
		x5d=x1*x1;
		if (step<2)
			goto L66;
		else
			goto L55;
L999:
		fclose (fptr);
		return 0;
}
