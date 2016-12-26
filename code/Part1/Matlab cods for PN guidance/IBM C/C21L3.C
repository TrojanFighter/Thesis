#include <stdio.h>
#include <stdlib.h>
#include <math.h>
main()
{
	int step;
	float vc,xnt,xnp,xnclim,x,w,tau,xmweaveold,xmweavemax;
	float tf,phase,y,yd,xnl,d,elamdh,x4,x5,t,h,s;
	float yold,ydold,xnlold,dold,elamdhold,x4old,x5old;
	float ydd,xnld,elamdhd,dd,x4d,x5d,ytdd,tgo,xlam,xnc;
	float xmweave;
	FILE *fptr;
	vc=4000.;
	xnt=32.2;
	xnp=3.;
	xnclim=99999.;
	printf("current values for target acceleration and acc limit %f\t%f\n",xnt,xnclim);
	printf ("enter target acceleration and acc limit:");
	scanf ("%f %f",&xnt,&xnclim);
	printf("you entered %f\t%f\n",xnt,xnclim);
	
	printf("current effective navigation ratio %f\n",xnp);
	printf ("enter effective navigation ratio:");
	scanf ("%f",&xnp);
	printf("you entered %f\n",xnp);
	
	fptr=fopen("DATFIL.TXT","w");
	for(x=.1;x<=4.;x=x+.1){
		if(x<.5){
			w=1.;
			tau=x/w;
		}
		else{
			w=x;
			tau=1.;
		}
		xmweaveold=0.;
		xmweavemax=0.;
		for(tf=.2;tf<=20.;tf=tf+.2){
			phase=0.;
			y=0.;
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
			s=s+h;
			goto L10;
L200:
			ytdd=xnt*sin(w*t);
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
			ydd=ytdd-xnl;
			if(step<=1)
				goto L66;
			else
				goto L55;
L999:
 			xmweave=y;
			if(xmweave>xmweaveold&&xmweave>xmweavemax&&tf>10.)
				xmweavemax=xmweave;
			xmweaveold=xmweave;
		}
 		if(x<.5)
			xmweavemax=xmweavemax/(tau*tau);
		printf("%10.3f %10.3f \n",x,xmweavemax);
		fprintf(fptr,"%10.3f %10.3f \n",x,xmweavemax);
	}
 	fclose (fptr);
	return 0;
}
