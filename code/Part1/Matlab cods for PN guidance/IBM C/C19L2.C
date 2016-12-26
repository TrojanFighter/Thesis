#include <stdio.h>
#include <stdlib.h>
#include <math.h>
main()
{
	float vc,xnt,xnclimg,yic,vm,hedeg,tau,xnp,ta,r,tf;
	float y,yd,ydic,xnl,elamdh,x4,x5,th,thh,t,h,s,xnclim;
	float yold,ydold,xnlold,elamdhold,x4old,x5old,thold;
	float thhold,ydd,x4d,x5d,thd,thhd,xnld,tgo,xlam,eps;
	float dd,xnc,elamdhd;
	int step;
	FILE *fptr;
	fptr=fopen("DATFIL.TXT","w");
	vc=4000.;
	xnt=32.2;
	xnclimg=7.;
	yic=0.;
	vm=3000.;
	hedeg=0.;
	tau=.3;
	xnp=3.;
	ta=5.;
	r=-.01;
	tf=10.;
	y=yic;
	yd=-vm*hedeg/57.3;
	ydic=yd;
	xnl=0.;
	elamdh=0.;
	x4=0.;
	x5=0.;
	th=0.;
	thh=0.;
	t=0.;
	h=.01;
	s=0.;
	xnclim=xnclimg*32.2;
L10:	if(t>(tf-.0001))goto L999;
 	yold=y;
	ydold=yd;
	xnlold=xnl;
	elamdhold=elamdh;
	x4old=x4;
	x5old=x5;
	thold=th;
	thhold=thh;
	step=1;
	goto L200;
L66:	step=2;
 	y=y+h*yd;
 	yd=yd+h*ydd;
	xnl=xnl+h*xnld;
	elamdh=elamdh+h*elamdhd;
	x4=x4+h*x4d;
	x5=x5+h*x5d;
	th=th+h*thd;
	thh=thh+h*thhd;
	t=t+h;
	goto L200;
L55:
 	y=.5*(yold+y+h*yd);
 	yd=.5*(ydold+yd+h*ydd);
	xnl=.5*(xnlold+xnl+h*xnld);
	elamdh=.5*(elamdhold+elamdh+h*elamdhd);
	x4=.5*(x4old+x4+h*x4d);
	x5=.5*(x5old+x5+h*x5d);
	th=.5*(thold+th+h*thd);
	thh=.5*(thhold+thh+h*thhd);
	s=s+h;
	if(s>.09999){
		s=0.;
		printf("%6.3f %6.3f %6.3f\n",t,y,xnc/32.2);
		fprintf(fptr,"%6.3f %6.3f %6.3f\n",t,y,xnc/32.2);
	}
	goto L10;
L200:
 	tgo=tf-t+.00001;
	xlam=y/(vc*tgo);
	eps=xlam-th-thh+r*thh;
	dd=5.*eps/tau;
	elamdhd=5.*(dd-elamdh)/tau;
	xnc=xnp*vc*elamdh;
	if(xnc>xnclim)xnc=xnclim;
	if(xnc<-xnclim)xnc=-xnclim;
	x4d=5.*(xnc-x4)/tau;
	x5d=5.*(x4-x5)/tau;
	xnld=5.*(x5-xnl)/tau;
	thd=xnl/vm+ta*xnld/vm;
	thhd=dd-thd;
	ydd=xnt-xnl;
	if (step<2)
		goto L66;
	else
		goto L55;
L999:	
	fclose (fptr);
	return 0;
}
	
