#include <stdio.h>
#include <stdlib.h>
#include <math.h>
main()
{
	int step;
	float xnt,xnp,tau,tf,vc,w,t,s,tp,x2,x3,x4,x5,x6;
	float h,x2old,x3old,x4old,x5old,x6old,x2d,x3d;
	float x4d,x5d,x6d,xmweave,y1,tgo;
	FILE *fptr;
	fptr=fopen("DATFIL.TXT","w");
	xnt=193.2;
	xnp=3.;
	tau=1.;
	tf=10.;
	vc=4000.;
	w=3.;
	t=0.;
	s=0.;
	tp=t+.00001;
	x2=0;
	x3=1;
	x4=0.;
	x5=0.;
	x6=0.;
	h=.01;
L10:
	if(tp>(tf-.00001))goto L999;
	s=s+h;
	x2old=x2;
	x3old=x3;
	x4old=x4;
	x5old=x5;
	x6old=x6;
	step=1;
	goto L200;
L66:
	step=2;
	x2=x2+h*x2d;
	x3=x3+h*x3d;
	x4=x4+h*x4d;
	x5=x5+h*x5d;
	x6=x6+h*x6d;
	tp=tp+h;
	goto L200;
L55:
	x2=(x2old+x2)/2+.5*h*x2d;
	x3=(x3old+x3)/2+.5*h*x3d;
	x4=(x4old+x4)/2+.5*h*x4d;
	x5=(x5old+x5)/2+.5*h*x5d;
	x6=(x6old+x6)/2+.5*h*x6d;
	if(s<.09999)goto L10;
	s=0.;
	xmweave=xnt*w*x6;
	printf("%10.3f %10.3f \n",tp,xmweave);
	fprintf(fptr,"%10.3f %10.3f \n",tp,xmweave);
	goto L10;
L200:
	x2d=x3;
	y1=(-x2+x4)/tau;
	tgo=tp+.00001;
	x3d=y1*xnp/tgo;
	x4d=-y1;
	x5d=x2-w*w*x6;
	x6d=x5;
	if(step<=1)
		goto L66;
	else
		goto L55;
L999:
	fclose (fptr);
	return 0;
}
