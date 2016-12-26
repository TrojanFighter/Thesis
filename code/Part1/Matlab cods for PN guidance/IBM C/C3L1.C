#include <stdio.h>
#include <stdlib.h>
#include <math.h>
main()
{
	int step;
	float xnt,xnp,tau,tf,vm,hedeg,t,s,tp,x1,x2,x3,x4,h,he;
	float x1old,x2old,x3old,x4old,x1d,x2d,x3d,x4d,xmnt,xmhe;
	float y1,tgo;
	FILE *fptr;
	fptr=fopen("DATFIL.TXT","w");
	xnt=96.6;
	xnp=4.;
	tau=1.;
	tf=10.;
	vm=3000.;
	hedeg=-20.;
	t=0.;
	s=0.;
	tp=t+.00001;
	x1=0;
	x2=0;
	x3=1;
	x4=0;
	h=.01;
	he=hedeg/57.3;
L10:
	if(tp>(tf-.00001))goto L999;
	s=s+h;
	x1old=x1;
	x2old=x2;
	x3old=x3;
	x4old=x4;
	step=1;
	goto L200;
L66:
	step=2;
	x1=x1+h*x1d;
	x2=x2+h*x2d;
	x3=x3+h*x3d;
	x4=x4+h*x4d;
	tp=tp+h;
	goto L200;
L55:
	x1=(x1old+x1)/2+.5*h*x1d;
	x2=(x2old+x2)/2+.5*h*x2d;
	x3=(x3old+x3)/2+.5*h*x3d;
	x4=(x4old+x4)/2+.5*h*x4d;
	if(s<.09999)goto L10;
	s=0.;
	xmnt=xnt*x1;
	xmhe=-vm*he*x2;
	printf("%10.4f %10.4f %10.4f \n",tp,xmnt,xmhe);
	fprintf(fptr,"%10.4f %10.4f %10.4f \n",tp,xmnt,xmhe);
	goto L10;
L200:
	x1d=x2;
	x2d=x3;
	y1=(x4-x2)/tau;
	tgo=tp+.00001;
	x3d=xnp*y1/tgo;
	x4d=-y1;
	if(step<=1)
		goto L66;
	else
		goto L55;
L999:
	fclose (fptr);
	return 0;
}
