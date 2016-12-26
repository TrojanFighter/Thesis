#include <stdio.h>
#include <stdlib.h>
#include <math.h>
main()
{
	int step;
	float xnt,xnp,tau,tf,vc,xntd,xntdd,t,s,tp;
	float x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,h;
	float x1old,x2old,x3old,x4old,x5old,x6old,x7old,x8old,x9old,x10old;
	float x1d,x2d,x3d,x4d,x5d,x6d,x7d,x8d,x9d,x10d;
	float xmnt,xmntd,xmntdd,y1,tgo;
	FILE *fptr;
	fptr=fopen("DATFIL.TXT","w");
	xnt=32.2;
	xnp=3.;
	tau=1.;
	tf=10.;
	vc=4000.;
	xntd=32.2;
	xntdd=32.2;
	t=0.;
	s=0.;
	tp=t+.00001;
	x1=0;
	x2=0;
	x3=1;
	x4=0;
	x5=0.;
	x6=0.;
	x7=0.;
	x8=0.;
	x9=0.;
	x10=0.;
	h=.01;
L10:
	if(tp>(tf-.00001))goto L999;
	s=s+h;
	x1old=x1;
	x2old=x2;
	x3old=x3;
	x4old=x4;
	x5old=x5;
	x6old=x6;
	x7old=x7;
	x8old=x8;
	x9old=x9;
	x10old=x10;
	step=1;
	goto L200;
L66:
	step=2;
	x1=x1+h*x1d;
	x2=x2+h*x2d;
	x3=x3+h*x3d;
	x4=x4+h*x4d;
	x5=x5+h*x5d;
	x6=x6+h*x6d;
	x7=x7+h*x7d;
	x8=x8+h*x8d;
	x9=x9+h*x9d;
	x10=x10+h*x10d;
	tp=tp+h;
	goto L200;
L55:
	x1=(x1old+x1)/2+.5*h*x1d;
	x2=(x2old+x2)/2+.5*h*x2d;
	x3=(x3old+x3)/2+.5*h*x3d;
	x4=(x4old+x4)/2+.5*h*x4d;
	x5=(x5old+x5)/2+.5*h*x5d;
	x6=(x6old+x6)/2+.5*h*x6d;
	x7=(x7old+x7)/2+.5*h*x7d;
	x8=(x8old+x8)/2+.5*h*x8d;
	x9=(x9old+x9)/2+.5*h*x9d;
	x10=(x10old+x10)/2+.5*h*x10d;
	if(s<.09999)goto L10;
	s=0.;
	xmnt=fabs(xnt*x1);
	xmntd=fabs(xntd*x9);
	xmntdd=fabs(xntdd*x10);
	printf("%10.4f %10.4f %10.4f %10.4f\n",tp,xmnt,xmntd,xmntdd);
	fprintf(fptr,"%10.4f %10.4f %10.4f %10.4f\n",tp,xmnt,xmntd,xmntdd);
	goto L10;
L200:
	x1d=x2;
	x2d=x3;
	y1=5.*(5.*x5/tau+x4)/tau;
	tgo=tp+.00001;
	x3d=y1/(vc*tgo);
	x4d=-y1;
	x5d=-5.*x5/tau+5.*x6*xnp*vc/tau;
	x6d=-5.*x6/tau+5.*x7/tau;
	x7d=-5.*x7/tau+5.*x8/tau;
	x8d=-5.*x8/tau-x2;
	x9d=x1;
	x10d=x9;
	if(step<=1)
		goto L66;
	else
		goto L55;
L999:
	fclose (fptr);
	return 0;
}
