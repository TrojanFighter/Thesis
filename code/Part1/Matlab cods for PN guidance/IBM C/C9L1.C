#include <stdio.h>
#include <stdlib.h>
#include <math.h>
main()
{
	float m11,m12,m13,m22,m23,m33,k1,k2,k3,vc,xnt,vm,hedeg,sigrin;
	float ts,tf,ts2,ts3,ts4,ts5,phin,rtm,signoise,sigpos,sign2;
	float p11,p12,p13,p22,p23,p33,t,h,s,tgo,siggl;
	FILE *fptr;
	fptr=fopen("DATFIL.TXT","w");
	siggl=0.;
	vc=4000.;
	xnt=96.6;
	vm=3000.;
	hedeg=20.;
	sigrin=.001;
	ts=.1;
	tf=10.;
	ts2=ts*ts;
	ts3=ts2*ts;
	ts4=ts3*ts;
	ts5=ts4*ts;
	phin=xnt*xnt/tf;
	rtm=vc*tf;
	signoise=sqrt(sigrin*sigrin+(siggl/rtm)*(siggl/rtm));
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
 L10:
 	if(t>(tf-.0001))goto L999;
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
	printf("%8.5f %8.5f %8.5f %8.5f\n",t,k1,k2,k3);
	fprintf(fptr,"%8.5f %8.5f %8.5f %8.5f\n",t,k1,k2,k3);
	t=t+ts;
	goto L10;
L999:
 	printf("%8.5f %8.5f %8.5f %8.5f\n",t,k1,k2,k3);
	fclose (fptr);
	return 0;
}
