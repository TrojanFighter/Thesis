#include <stdio.h>
#include <stdlib.h>
#include <math.h>
void distance(float x,float y,float xfirst,float yfirst,float *distnm);
main()
{
	int step,left;
	float xisp1,xisp2,xmf1,xmf2,wpay,delv,delv1,delv2,amax1,amax2;
	float gamdeg,top2,bot2,wp2,ws2,wtot2,trst2,tb2,top1,bot1;
	float wp1,ws1,wtot,trst1,tb1,delvk;
	float h,t,s,a,gm,altnm,angdeg,ang,vrx,vry,x,y,alt,xfirst,yfirst;
	float x1,y1,xold,yold,x1old,y1old,xd,yd,x1d,y1d,distnm,xnm,ynm;
	float wgt,trst,at,vel,axt,ayt,tembot;
	FILE *fptr;
	fptr=fopen("DATFIL.TXT","w");
	left=1;
	xisp1=250.;
	xisp2=250.;
	xmf1=.85;
	xmf2=.85;
	wpay=100.;
	delv=20000.;
	delv1=.3333*delv;
	delv2=.6667*delv;
	amax1=10.;
	amax2=10.;
	gamdeg=85.;
	top2=wpay*(exp(delv2/(xisp2*32.2))-1.);
	bot2=1/xmf2-((1.-xmf2)/xmf2)*exp(delv2/(xisp2*32.2));
	wp2=top2/bot2;
	ws2=wp2*(1-xmf2)/xmf2;
	wtot2=wp2+ws2+wpay;
	trst2=amax2*(wpay+ws2);
	tb2=xisp2*wp2/trst2;
	top1=wtot2*(exp(delv1/(xisp1*32.2))-1.);
	bot1=1/xmf1-((1.-xmf1)/xmf1)*exp(delv1/(xisp1*32.2));
	wp1=top1/bot1;
	ws1=wp1*(1-xmf1)/xmf1;
	wtot=wp1+ws1+wtot2;
	trst1=amax1*(wtot2+ws1);
	tb1=xisp1*wp1/trst1;
	delvk=delv/1000.;
	h=.01;
	t=0.;
	s=0.;
	a=2.0926e7;
	gm=1.4077e16;
	altnm=0.;
	alt=altnm*6076.;
	angdeg=90.;
	ang=angdeg/57.3;
	if(left==1){
		vrx=cos(1.5708-gamdeg/57.3+ang);
		vry=sin(1.5708-gamdeg/57.3+ang);
	}
	else{
		vrx=cos(-1.5708+gamdeg/57.3+ang);
	  	vry=sin(-1.5708+gamdeg/57.3+ang);
	}
	x=(a+alt)*cos(ang);
	y=(a+alt)*sin(ang);
	alt=sqrt(x*x+y*y)-a;
	xfirst=x;
	yfirst=y;
	x1=vrx;
	y1=vry;
L10:	if(alt<0.&&t>10.)goto L999;
	xold=x;
	yold=y;
	x1old=x1;
	y1old=y1;
	step=1;
	goto L200;
L66:	step=2;
	x=x+h*xd;
	y=y+h*yd;
	x1=x1+h*x1d;
	y1=y1+h*y1d;
	t=t+h;
	goto L200;
L55:
	x=(xold+x)/2+.5*h*xd;
	y=(yold+y)/2+.5*h*yd;
	x1=(x1old+x1)/2+.5*h*x1d;
	y1=(y1old+y1)/2+.5*h*y1d;
	alt=sqrt(x*x+y*y)-a;
 	s=s+h;
	if(s<9.99999)goto L10;
 	s=0.;
 	distance(x,y,xfirst,yfirst,&distnm);
	altnm=(sqrt(x*x+y*y)-a)/6076.;
	xnm=x/6076.;
	ynm=y/6076.;
	printf("%10.3f %10.3f %10.3f %10.3f %10.3f\n",t,distnm,altnm,xnm,ynm);
	fprintf(fptr,"%10.3f %10.3f %10.3f %10.3f %10.3f\n",t,distnm,altnm,xnm,ynm);
	goto L10;
L200:
 	if(t<tb1){
		wgt=-wp1*t/tb1+wtot;
		trst=trst1;
	}
	else if(t<(tb1+tb2)){
		wgt=-wp2*t/tb2+wtot2+wp2*tb1/tb2;
		trst=trst2;
	}
	else{
		wgt=wpay;
		trst=0.;
	}
	at=32.2*trst/wgt;
	vel=sqrt(x1*x1+y1*y1);
	axt=at*x1/vel;
	ayt=at*y1/vel;
	tembot=pow(x*x+y*y,1.5);
	x1d=-gm*x/tembot+axt;
	y1d=-gm*y/tembot+ayt;
	xd=x1;
	yd=y1;
	if(step<=1)
		goto L66;
	else
		goto L55;
L999:
 	distance(x,y,xfirst,yfirst,&distnm);
	altnm=(sqrt(x*x+y*y)-a)/6076.;
	xnm=x/6076.;
	ynm=y/6076.;
	printf("%10.3f %10.3f %10.3f %10.3f %10.3f\n",t,distnm,altnm,xnm,ynm);
	fprintf(fptr,"%10.3f %10.3f %10.3f %10.3f %10.3f\n",t,distnm,altnm,xnm,ynm);
	fclose (fptr);
	return 0;
}
	
	void distance(float x,float y,float xfirst,float yfirst,float *distnm)
{
	float r,a,cbeta,beta;
	r=sqrt(x*x+y*y);
	a=2.0926E7;
	cbeta=(x*xfirst+y*yfirst)/(r*a);
	if(cbeta<1.){
		beta=atan2(sqrt(1.-cbeta*cbeta),cbeta);
		*distnm=a*beta/6076.;
	}
	else
		*distnm=(xfirst-x)/6076.;
}
