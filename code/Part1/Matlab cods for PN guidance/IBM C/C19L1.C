#include <stdio.h>
#include <stdlib.h>
#include <math.h>
main()
{
	int step;
	float xnp,tau,tf,vc,phifn,phirn,phirna,phigl,ra;
	float t,s,tp,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,h;
	float x2old,x3old,x4old,x5old,x6old,x7old,x8old;
	float x9old,x10old,x11old,x12old,temp1;
	float x2d,x3d,x4d,x5d,x6d,x7d,x8d,x9d,x10d,x11d,x12d;
	float xmfn,xmrn,xmrna,xmgl,y1,tgo,phin;
	FILE *fptr;
	fptr=fopen("DATFIL.TXT","w");
	xnp=3.;
	tau=1.;
	tf=10.;
	vc=1.;
	phin=1.;
	phirn=1.;
	phirna=1.;
	phigl=1.;
	phifn=1.;
	ra=1.;
	t=0.;
	s=0.;
	tp=t+.00001;
	x2=0;
	x3=1;
	x4=0;
	x5=0.;
	x6=0.;
	x7=0.;
	x8=0.;
	x9=0.;
	x10=0.;
	x11=0.;
	x12=0.;
	h=.01;
L10:	if(tp>(tf-.00001))goto L999;
	s=s+h;
	x2old=x2;
	x3old=x3;
	x4old=x4;
	x5old=x5;
	x6old=x6;
	x7old=x7;
	x8old=x8;
	x9old=x9;
	x10old=x10;
	x11old=x11;
	x12old=x12;
	step=1;
	goto L200;
L66:	step=2;
	x2=x2+h*x2d;
	x3=x3+h*x3d;
	x4=x4+h*x4d;
	x5=x5+h*x5d;
	x6=x6+h*x6d;
	x7=x7+h*x7d;
	x8=x8+h*x8d;
	x9=x9+h*x9d;
	x10=x10+h*x10d;
	x11=x11+h*x11d;
	x12=x12+h*x12d;
	tp=tp+h;
	goto L200;
L55:
	x2=(x2old+x2)/2+.5*h*x2d;
	x3=(x3old+x3)/2+.5*h*x3d;
	x4=(x4old+x4)/2+.5*h*x4d;
	x5=(x5old+x5)/2+.5*h*x5d;
	x6=(x6old+x6)/2+.5*h*x6d;
	x7=(x7old+x7)/2+.5*h*x7d;
	x8=(x8old+x8)/2+.5*h*x8d;
	x9=(x9old+x9)/2+.5*h*x9d;
	x10=(x10old+x10)/2+.5*h*x10d;
	x11=(x11old+x11)/2+.5*h*x11d;
	x12=(x12old+x12)/2+.5*h*x12d;
	if(s<.09999)goto L10;
	s=0.;
	xmfn=sqrt(x9*phifn);
	xmrn=sqrt(x10*phirn);
	xmrna=sqrt(x11*phirna);
	xmgl=sqrt(x12*phigl);
	printf("%6.3f %6.3f %6.3f %6.3f %6.3f\n",tp,xmfn,xmrn,xmrna,xmgl);
	fprintf(fptr,"%6.3f %6.3f %6.3f %6.3f %6.3f\n",tp,xmfn,xmrn,xmrna,xmgl);
	goto L10;
L200:
	x2d=x3;
	y1=5.*(5.*x5/tau+x4)/tau;
	tgo=tp+.00001;
	x3d=y1/(vc*tgo);
	x4d=-y1;
	x5d=-5.*x5/tau+5.*x6*xnp*vc/tau;
	x6d=-5.*x6/tau+5.*x7/tau;
	x7d=-5.*x7/tau+5.*x8/tau;
	x8d=-5.*x8/tau-x2;
	x9d=y1*y1;
	x10d=(y1*vc*tgo/ra)*(y1*vc*tgo/ra);
	temp1=(vc*tgo/ra);
	x11d=(y1*(temp1*temp1))*(y1*(temp1*temp1));
	x12d=(y1/(vc*tgo))*(y1/(vc*tgo));
	if (step<2)
		goto L66;
	else
		goto L55;
L999:
	fclose (fptr);
	return 0;
}
