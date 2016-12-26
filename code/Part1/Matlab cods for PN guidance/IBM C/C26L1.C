#include <stdio.h>
#include <stdlib.h>
#include <math.h>
void gauss(float signoise,float *xlamnoise);
main()
{	
	float m11,m12,m13,m22,m23,m33,k1,k2,k3,tau,vc,xnt,xntreal;
	float xntmax,w,yic,vm,hedeg,hedegfil,xnp,ts,tf,y,yd,ydic;
	float ts2,ts3,ts4,ts5,phin,rtm,signoise,sigpos,sign2,p11;
	float p12,p13,p22,p23,p33,t,h,s,yh,ydh,xnth,xnc,xnl,yold;
	float ydold,xnlold,xlamd,xnld,ydd,sigrin,xlamnoise;
	float tgo,ystar,res,xlamdh,errnt,sp33,sp33p,theory,xlam;
	float xs,top,bot1,bot2,xnpp,c1,c2,c3,c4;
	int step,apn;
	FILE *fptr;
	fptr=fopen("DATFIL.TXT","w");
	tau=.5;
	apn=0;
	vc=9000.;
	xnt=96.6;
	xntreal=96.6;
	xntmax=96.6;
	w=2.;
	yic=0.;
	vm=3000.;
	hedeg=0.;
	hedegfil=20.;
	xnp=3.;
	sigrin=.001;
	ts=.01;
	tf=10.;
	y=yic;
	yd=-vm*hedeg/57.3;
	ydic=yd;
	ts2=ts*ts;
	ts3=ts2*ts;
	ts4=ts3*ts;
	ts5=ts4*ts;
	phin=xntmax*xntmax/tf;
	rtm=vc*tf;
	signoise=sigrin;
	sigpos=rtm*signoise;
	sign2=sigpos*sigpos;
	p11=sign2;
	p12=0.;
	p13=0.;
	p22=(vm*hedegfil/57.3)*(vm*hedegfil/57.3);
	p23=0.;
	p33=xntmax*xntmax;
	t=0.;
	h=.001;
	s=0.;
	yh=0.;
	ydh=0.;
	xnth=0.;
	xnc=0.;
	xnl=0.;
L10:
	if(t>(tf-.0001))goto L999;
 	yold=y;
	ydold=yd;
	xnlold=xnl;
	step=1;
	goto L200;
L66:
	step=2;
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
	signoise=sigrin;
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
	gauss(signoise,&xlamnoise);
	ystar=rtm*(xlam+xlamnoise);
	res=ystar-yh-ts*ydh-.5*ts*ts*(xnth-xnl);
	yh=k1*res+yh+ts*ydh+.5*ts*ts*(xnth-xnl);
	ydh=k2*res+ydh+ts*(xnth-xnl);
	xnth=k3*res+xnth;
	xlamdh=(yh+ydh*tgo)/(vc*tgo*tgo);
	if(apn==0){
		xnc=xnp*(yh+ydh*tgo)/(tgo*tgo);
	}
	if(apn==1){
		xnc=xnp*(yh+ydh*tgo+.5*xnth*tgo*tgo)/(tgo*tgo);
	}
	if(apn==2){
		xs=tgo/tau;
		top=6.*xs*xs*(exp(-xs)-1.+xs);
		bot1=2*xs*xs*xs+3.+6.*xs-6.*xs*xs;
		bot2=-12.*xs*exp(-xs)-3.*exp(-2.*xs);
		xnpp=top/(.0001+bot1+bot2);
		c1=xnpp/(tgo*tgo);
		c2=xnpp/tgo;
		c3=.5*xnpp;
		c4=-xnpp*(exp(-xs)+xs-1.)/(xs*xs);
		xnc=c1*yh+c2*ydh+c3*xnth+c4*xnl;
	}
	errnt=xnt-xnth;
	sp33=sqrt(p33);
	sp33p=-sp33;
	theory=sqrt(m11+sign2);
	printf("%6.3f %6.3f %6.3f %6.3f %6.3f\n",t,xnt/32.2,xnth/32.2,y,ystar);
	fprintf(fptr,"%6.3f %6.3f %6.3f %6.3f %6.3f\n",t,xnt/32.2,xnth/32.2,y,ystar);
	goto L10;
L200:
	xnt=xntreal*sin(w*t);
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
	printf("%6.3f\n",y);
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
