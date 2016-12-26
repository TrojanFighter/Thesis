#include <stdio.h>
#include <stdlib.h>
#include <math.h>
void uniform(float *sum);
main()
{
	float z[1000],vc,xnt,vm,xnp,tau,tf,z1,sum,tstart,pz;
	float coef,y,yd,t,h,s,xnc,xnl;
	float yold,ydold,xnlold,ydd,xnld,tgo,rtm,xlamd,xmean,sigma;
	int run,i,step;
	FILE *fptr;
	fptr=fopen("DATFIL.TXT","w");
	vc=4000.;
	xnt=96.6;
	vm=3000.;
	xnp=3.;
	tau=1.;
	run=100;
L106:
	for (tf=1.; tf<=10.; tf=tf+.5){
		z1=0.;
		for (i=1; i<=run; i=i+1){
			uniform(&sum);
			tstart=tf*sum;
			uniform(&pz);
			pz=pz-.5;
			if(pz>0.)
				coef=1.;
			else
				coef=-1.;
			y=0.;
			yd=0.;
			t=0.;
			h=.01;
			s=0.;
			xnc=0.;
			xnl=0.;
L10:			if(t>(tf-.0001))goto L999;
 			if(t<tstart)
				xnt=0.;
			else
				xnt=coef*96.6;
 			yold=y;
			ydold=yd;
			xnlold=xnl;
			step=1;
			goto L200;
L66:			step=2;
 			y=y+h*yd;
 			yd=yd+h*ydd;
			xnl=xnl+h*xnld;
			t=t+h;
			goto L200;
L55:
 			y=.5*(yold+y+h*yd);
 			yd=.5*(ydold+yd+h*ydd);
			xnl=.5*(xnlold+xnl+h*xnld);
			s=s+h;
			goto L10;
L200:
 			tgo=tf-t+.00001;
			rtm=vc*tgo;
			xlamd=(rtm*yd+y*vc)/(rtm*rtm);
			xnc=xnp*vc*xlamd;
			xnld=(xnc-xnl)/tau;
			ydd=xnt-xnl;
			if (step<2)
				goto L66;
			else
				goto L55;
L999:
			z[i]=y;
			z1=z[i]+z1;
			xmean=z1/i;
		}
 		sigma=0.;
		z1=0.;
		for (i=1; i<=run; i=i+1){
			z1=(z[i]-xmean)*(z[i]-xmean)+z1;
			if(i==1)
				sigma=0.;
			else
				sigma=sqrt(z1/(i-1));
 		}
 		printf("%6.3f %6.3f %6.3f\n",tf,sigma,xmean);
		fprintf(fptr,"%6.3f %6.3f %6.3f\n",tf,sigma,xmean);
 	}
 	fclose (fptr);
 	return 0;
}

	void uniform(float *sum)
{
	int x1;
	x1=rand();
	*sum=(float)x1/(float)RAND_MAX;
}
