#include <stdio.h>
#include <stdlib.h>
#include <math.h>
main()
{
	int step;
	float vc,xnt,displace,vm,tau,xnp,xnclim,y,yd,xnl,d,elamdh,x4;
	float x5,t,h,s,yold,ydold,xnlold,dold,elamdhold,x4old,x5old;
	float ydd,xnld,elamdhd,dd,x4d,x5d,tgo,xlam,xnc,tf;
	FILE *fptr;
	fptr=fopen("DATFIL.TXT","w");
	vc=4000.;
	xnt=0.;
	displace=1.;
	vm=3000.;
	tau=1.;
	xnp=3.;
	xnclim=99999999.;
	for(tf=.1;tf<=10.;tf=tf+.1){
		y=displace;
		yd=0.;
		xnl=0.;
		d=0.;
		elamdh=0.;
		x4=0.;
		x5=0.;
		t=0.;
		h=.01;
		s=0.;
L10:
		if(t>(tf-.0001))goto L999;
 		yold=y;
		ydold=yd;
		xnlold=xnl;
		dold=d;
		elamdhold=elamdh;
		x4old=x4;
		x5old=x5;
		step=1;
		goto L200;
L66:
		step=2;
 		y=y+h*yd;
 		yd=yd+h*ydd;
		xnl=xnl+h*xnld;
		elamdh=elamdh+h*elamdhd;
		d=d+h*dd;
		x4=x4+h*x4d;
		x5=x5+h*x5d;
		t=t+h;
		goto L200;
L55:
 		y=.5*(yold+y+h*yd);
 		yd=.5*(ydold+yd+h*ydd);
		xnl=.5*(xnlold+xnl+h*xnld);
		d=.5*(dold+d+h*dd);
		elamdh=.5*(elamdhold+elamdh+h*elamdhd);
		x4=.5*(x4old+x4+h*x4d);
		x5=.5*(x5old+x5+h*x5d);
		goto L10;
L200:
		tgo=tf-t+.00001;
		xlam=y/(vc*tgo);
		dd=5.*(xlam-d)/tau;
		elamdhd=5.*(dd-elamdh)/tau;
		xnc=xnp*vc*elamdh;
		if(xnc>xnclim)xnc=xnclim;
		if(xnc<-xnclim)xnc=-xnclim;
		x4d=5.*(xnc-x4)/tau;
		x5d=5.*(x4-x5)/tau;
		xnld=5.*(x5-xnl)/tau;
		ydd=xnt-xnl;
		if(step<=1)
			goto L66;
		else
			goto L55;
L999:
		printf("%10.3f %10.3f \n",tf,y);
		fprintf(fptr,"%10.3f %10.3f \n",tf,y);
 	}
	fclose (fptr);
	return 0;
}
