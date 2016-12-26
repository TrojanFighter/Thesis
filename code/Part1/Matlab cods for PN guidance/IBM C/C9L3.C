#include <stdio.h>
#include <stdlib.h>
#include <math.h>
void gauss(float signoise,float *xlamnoise);
main()
{
	float z[1000],vc,xnt,yic,vm,hedeg,xnp,signoise;
	int run,apn,noise,step,i;
	float m11,m12,m13,m22,m23,m33,k1,k2,k3;
	float ts,tau,xlim,tf,z1,y,yd,ts2,ts3,ts4,ts5;
	float phin,rtm,sigpos,sign2,p11,p12,p13,p22,p23,p33;
	float t,h,s,yh,ydh,xnth,xnc,xnl,yold,ydold,xnlold;
	float ydd,xnld,tgo,ystar,res,xlamdh,x,top,bot1,bot2,xnpp,xnew;
	float xlam,xlamd,xmean,sigma,xlamnoise;
	FILE *fptr;
	fptr=fopen("DATFIL.TXT","w");
	vc=4000.;
	xnt=96.6;
	yic=0.;
	vm=3000.;
	hedeg=20.;
	xnp=3.;
	signoise=.001;
	ts=.1;
	tau=.5;
	noise=1;
	run=50;
	apn=0;
	xlim=999999.;
 L106:
	for (tf=.5; tf<=10.; tf=tf+.5){
		z1=0.;
		for (i=1; i<=run; i=i+1){
			y=yic;
			yd=0.;
			ts2=ts*ts;
			ts3=ts2*ts;
			ts4=ts3*ts;
			ts5=ts4*ts;
			phin=xnt*xnt/tf;
			rtm=vc*tf;
			sigpos=rtm*signoise;
			sign2=sigpos*sigpos;
			p11=sign2;
			p12=0.;
			p13=0.;
			p22=(vm*hedeg/57.3)*(vm*hedeg/57.3);
			p23=0.;
			p33=xnt*xnt;
			t=0.;
			h=.01;
			s=0.;
			yh=0.;
			ydh=0.;
			xnth=0.;
			xnc=0.;
			xnl=0.;
L10:			if(t>(tf-.0001))goto L999;
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
			if(s<(ts-.0001))goto L10;
			s=0.;
			tgo=tf-t+.000001;
			rtm=vc*tgo;
			sigpos=rtm*signoise;
			sign2=sigpos*sigpos;
			m11=p11+ts*p12+.5*ts2*p13+ts*(p12+ts*p22+.5*ts2*p23);
			m11=m11+.5*ts2*(p13+ts*p23+.5*ts2*p33)+ts5*phin/20.;
			m12=p12+ts*p22+.5*ts2*p23+ts*(p13+ts*p23+.5*ts2*p33)+ts4*phin/8.;
			m13=p13+ts*p23+.5*ts2*p33+phin*ts3/6.;
			m22=p22+ts*p23+ts*(p23+ts*p33)+phin*ts3/3.;
			m23=p23+ts*p33+.5*ts2*phin;
			m33=p33+phin*ts;
			k1=m11/(m11+sign2);
			k2=m12/(m11+sign2);
			k3=m13/(m11+sign2);
			p11=(1.-k1)*m11;
			p12=(1.-k1)*m12;
			p13=(1.-k1)*m13;
			p22=-k2*m12+m22;
			p23=-k2*m13+m23;
			p33=-k3*m13+m33;
			if(noise==1)
				gauss(signoise,&xlamnoise);
			else
				xlamnoise=0.;
			ystar=rtm*(xlam+xlamnoise);
			res=ystar-yh-ts*ydh-.5*ts*ts*(xnth-xnl);
			yh=k1*res+yh+ts*ydh+.5*ts*ts*(xnth-xnl);
			ydh=k2*res+ydh+ts*(xnth-xnl);
			xnth=k3*res+xnth;
			xlamdh=(yh+ydh*tgo)/(vc*tgo*tgo);
			if(apn==0)
				xnc=xnp*vc*xlamdh;
			else if(apn==1)
				xnc=xnp*vc*xlamdh+apn*.5*xnp*xnth;
			else{
				x=tgo/tau;
				top=6.*x*x*(exp(-x)-1.+x);
				bot1=2*x*x*x+3.+6.*x-6.*x*x;
				bot2=-12.*x*exp(-x)-3.*exp(-2.*x);
				xnpp=top/(.0001+bot1+bot2);
				xnew=xnpp*xnl*(exp(-x)+x-1.)/(x*x);
				xnc=xnpp*vc*xlamdh+.5*xnpp*xnth-xnew;
			}
			if(xnc>xlim)xnc=xlim;
			if(xnc<-xlim)xnc=-xlim;
			goto L10;
L200:	
 			tgo=tf-t+.00001;
			rtm=vc*tgo;
			xlam=y/(vc*tgo);
			xlamd=(rtm*yd+y*vc)/(rtm*rtm);
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
	
	void gauss(float signoise,float *xlamnoise)
{
	float sum,x,y,temp;
	int j;
	sum=0.;
	for (j=1; j<=6; j=j+1){
		x=rand();
		y=(float)x/(float)RAND_MAX;
		sum=sum+y;
		}
	temp=sum-3.;
	*xlamnoise=1.414*temp*signoise;
}
