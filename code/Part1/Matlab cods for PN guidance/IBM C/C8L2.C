#include <stdio.h>
#include <stdlib.h>
#include <math.h>
main()
{
	float xnt,xnp,tau,tf,vm,hedeg,t,s,tp,x1,x2,x3,x4,xnpp,h,he;
	float x1old,x2old,x3old,x4old,x1d,x2d,x3d,x4d,xmnt,xmhe;
	float c1,c2,c3,c4,x,top,bot1,bot2,tgo;
	int step,apn;
	FILE *fptr;
	fptr=fopen("DATFIL.TXT","w");
	xnt=96.6;
	xnp=4.;
	tau=1.;
	tf=10.;
	vm=3000.;
	hedeg=-20.;
	apn=0;
	t=0.;
	s=0.;
	tp=t+.00001;
	x1=0.;
	x2=0.;
	x3=1.;
	x4=0.;
	xnpp=0.;
	h=.01;
	he=hedeg/57.3;
L10:
	if(tp>(tf-.00001))
		goto L999;
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
	if(s<.09999)
		goto L10;
	s=0.;
	xmnt=xnt*x1;
	xmhe=-vm*he*x2;
	printf("%8.5f %8.5f %8.5f\n", tp,xmnt,xmhe);
	fprintf(fptr,"%8.5f %8.5f %8.5f\n", tp,xmnt,xmhe);
	goto L10;
L200:
 	tgo=tp+.00001;
	if(apn==0){
		c1=xnp/(tgo*tgo);
		c2=xnp/tgo;
		c3=0.;
		c4=0.;
	}
	if(apn==1){
		c1=xnp/(tgo*tgo);
		c2=xnp/tgo;
		c3=.5*xnp;
		c4=0.;
	}
	if(apn==2){
		x=tgo/tau;
		top=6.*x*x*(exp(-x)-1.+x);
		bot1=2*x*x*x+3.+6.*x-6.*x*x;
		bot2=-12.*x*exp(-x)-3.*exp(-2.*x);
		xnpp=top/(.0001+bot1+bot2);
		c1=xnpp/(tgo*tgo);
		c2=xnpp/tgo;
		c3=.5*xnpp;
		c4=-xnpp*(exp(-x)+x-1.)/(x*x);
	}
	x1d=x2+c3*x4/tau;
	x2d=x3+c2*x4/tau;
	x3d=c1*x4/tau;
	x4d=-x4/tau-x2+c4*x4/tau;
	if(step<=1)
		goto L66;
	else
		goto L55;
L999:
	fclose (fptr);
	return 0;
}
	
