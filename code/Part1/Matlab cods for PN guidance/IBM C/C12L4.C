#include <stdio.h>
#include <stdlib.h>
#include <math.h>
void distance(float x,float y,float xfirst,float yfirst,float *distnm);
main()
{
	int step;
	float h,a,gm,gamdeg,gam,distnm,angdeg,ang,phi,altnm,alt;
	float r0,top,temp,bot,v,vrx,vry,s,scount,x,y,xfirst,yfirst;
	float x1,y1,t,xold,yold,x1old,y1old,xd,yd,x1d,y1d,xnm,ynm,tembot;
	FILE *fptr;
	fptr=fopen("DATFIL.TXT","w");
 	h=.01;
	a=2.0926e7;
	gm=1.4077e16;
	gamdeg=23.;
	gam=gamdeg/57.3;
	distnm=6000.;
	angdeg=0.;
	ang=angdeg/57.3;
	phi=distnm*6076./a;
	altnm=0.;
	alt=altnm*6076.;
	r0=a+alt;
	top=gm*(1.-cos(phi));
	temp=r0*cos(gam)/a-cos(phi+gam);
	bot=r0*cos(gam)*temp;
	v=sqrt(top/bot);
	vrx=v*cos(1.5708-gam+ang);
	vry=v*sin(1.5708-gam+ang);
	s=0.;
	scount=0.;
	x=(a+alt)*cos(ang);
	y=(a+alt)*sin(ang);
	xfirst=x;
	yfirst=y;
	x1=vrx;
	y1=vry;
	t=0.;
L10:	if(alt<0.)goto L999;
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
L55:	x=(xold+x)/2+.5*h*xd;
	y=(yold+y)/2+.5*h*yd;
	x1=(x1old+x1)/2+.5*h*x1d;
	y1=(y1old+y1)/2+.5*h*y1d;
	alt=sqrt(x*x+y*y)-a;
	s=s+h;
	scount=scount+h;
	if(scount<9.99999)goto L10;
 	scount=0.;
	xnm=x/6076.;
	ynm=y/6076.;
	altnm=alt/6076.;
	distance(x,y,xfirst,yfirst,&distnm);
	printf("%8.2f %8.2f %8.2f %8.2f %8.2f\n",t,xnm,ynm,distnm,altnm);
	fprintf(fptr,"%8.2f %8.2f %8.2f %8.2f %8.2f\n",t,xnm,ynm,distnm,altnm);
 	goto L10;
L200:
	tembot=pow(x*x+y*y,1.5);
	x1d=-gm*x/tembot;
	y1d=-gm*y/tembot;
	xd=x1;
	yd=y1;
	if(step<=1)
		goto L66;
	else
		goto L55;
L999:
 	xnm=x/6076.;
	ynm=y/6076.;
	altnm=alt/6076.;
	distance(x,y,xfirst,yfirst,&distnm);
	printf("%8.2f %8.2f %8.2f %8.2f %8.2f\n",t,xnm,ynm,distnm,altnm);
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
