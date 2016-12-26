#include <stdio.h>
#include <stdlib.h>
#include <math.h>
void distance(float x,float y,float xfirst,float yfirst,float *distnm);
main()
{
	int step;
	float h,a,gm,gam,altnm,v,angdeg,ang,vrx,vry,alt,s,scount;
	float r0,r1,psi,x,y,xfirst,yfirst,x1,y1,t;
	float r0old,r1old,psiold,xold,yold,x1old,y1old;
	float r0d,r1d,psid,xd,yd,x1d,y1d,spolarnm,altpolarnm,distnm,tembot;
	FILE *fptr;
	fptr=fopen("DATFIL.TXT","w");
 	h=.01;
	a=2.0926e7;
	gm=1.4077e16;
	gam=45.;
	altnm=0.;
	v=24000.;
	angdeg=0.;
	ang=angdeg/57.3;
	vrx=v*cos(1.5708-gam/57.3+ang);
	vry=v*sin(1.5708-gam/57.3+ang);
	alt=altnm/6076.;
	s=0.;
	scount=0.;
	r0=a+alt;
	r1=v*sin(gam/57.3);
	psi=0.;
	x=(a+alt)*cos(ang);
	y=(a+alt)*sin(ang);
	xfirst=x;
	yfirst=y;
	x1=vrx;
	y1=vry;
	t=0.;
L10:	if(altnm<0.)goto L999;
 	r0old=r0;
	r1old=r1;
	psiold=psi;
	xold=x;
	yold=y;
	x1old=x1;
	y1old=y1;
	step=1;
	goto L200;
L66:	step=2;
 	r0=r0+h*r0d;
	r1=r1+h*r1d;
	psi=psi+h*psid;
	x=x+h*xd;
	y=y+h*yd;
	x1=x1+h*x1d;
	y1=y1+h*y1d;
	t=t+h;
	goto L200;
L55:	r0=(r0old+r0)/2+.5*h*r0d;
	r1=(r1old+r1)/2+.5*h*r1d;
	psi=(psiold+psi)/2+.5*h*psid;
	x=(xold+x)/2+.5*h*xd;
	y=(yold+y)/2+.5*h*yd;
	x1=(x1old+x1)/2+.5*h*x1d;
	y1=(y1old+y1)/2+.5*h*y1d;
	s=s+h;
	scount=scount+h;
	if(scount<9.99999)goto L10;
 	scount=0.;
	spolarnm=a*psi/6076.;
	altpolarnm=(r0-a)/6076.;
	altnm=(sqrt(x*x+y*y)-a)/6076.;
	distance(x,y,xfirst,yfirst,&distnm);
	printf("%8.2f %8.2f %8.2f %8.2f %8.2f\n",t,spolarnm,altpolarnm,distnm,altnm);
	fprintf(fptr,"%8.2f %8.2f %8.2f %8.2f %8.2f\n",t,spolarnm,altpolarnm,distnm,altnm);
 	goto L10;
 L200:
	psid=(a+alt)*v*cos(gam/57.3)/(r0*r0);
	r1d=-gm/(r0*r0)+r0*psid*psid;
	r0d=r1;
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
