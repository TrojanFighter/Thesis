#include <stdio.h>
#include <stdlib.h>
#include <math.h>
void distance(float x,float y,float xfirst,float yfirst,float *distnm);
main()
{
	int step;
	float h,a,gm,gam,altnm,v,alt,ang,vrx,vry,g,s,scount,rt1,rt2;
	float vt1,vt2,x,y,xfirst,yfirst,x1,y1,t;
	float rt1old,rt2old,vt1old,vt2old,xold,yold,x1old,y1old;
	float rt1d,rt2d,vt1d,vt2d,xd,yd,x1d,y1d;
	float rt1nm,rt2nm,distnm;
	float at1,at2,tembot;
	FILE *fptr;
	fptr=fopen("DATFIL.TXT","w");
 	h=.01;
	a=2.0926e7;
	gm=1.4077e16;
	gam=45.;
	altnm=0.;
	v=3000.;
	alt=altnm/6076.;
	ang=0.;
	vrx=v*cos(1.5708-gam/57.3+ang);
	vry=v*sin(1.5708-gam/57.3+ang);
	g=32.2;
	s=0.;
	scount=0.;
	rt1=alt*cos(ang);
	rt2=alt*sin(ang);
	vt1=vrx;
	vt2=vry;
	x=(a+alt)*cos(ang);
	y=(a+alt)*sin(ang);
	xfirst=x;
	yfirst=y;
	x1=vrx;
	y1=vry;
	t=0.;
L10:	if(altnm<0.)goto L999;
 	rt1old=rt1;
 	rt2old=rt2;
	vt1old=vt1;
	vt2old=vt2;
	xold=x;
	yold=y;
	x1old=x1;
	y1old=y1;
	step=1;
	goto L200;
L66:	step=2;
 	rt1=rt1+h*rt1d;
	rt2=rt2+h*rt2d;
	vt1=vt1+h*vt1d;
	vt2=vt2+h*vt2d;
	x=x+h*xd;
	y=y+h*yd;
	x1=x1+h*x1d;
	y1=y1+h*y1d;
	t=t+h;
	goto L200;
L55:	rt1=(rt1old+rt1)/2.+.5*h*rt1d;
 	rt2=(rt2old+rt2)/2.+.5*h*rt2d;
	vt1=(vt1old+vt1)/2.+.5*h*vt1d;
	vt2=(vt2old+vt2)/2.+.5*h*vt2d;
	x=(xold+x)/2+.5*h*xd;
	y=(yold+y)/2+.5*h*yd;
	x1=(x1old+x1)/2+.5*h*x1d;
	y1=(y1old+y1)/2+.5*h*y1d;
	s=s+h;
	scount=scount+h;
	if(scount<1.99999)goto L10;
 	scount=0.;
	rt1nm=rt1/6076.;
	rt2nm=rt2/6076.;
	altnm=(sqrt(x*x+y*y)-a)/6076.;
	distance(x,y,xfirst,yfirst,&distnm);
	printf("%8.2f %8.2f %8.2f %8.2f %8.2f\n",t,rt1nm,rt2nm,distnm,altnm);
	fprintf(fptr,"%8.2f %8.2f %8.2f %8.2f %8.2f\n",t,rt1nm,rt2nm,distnm,altnm);
 	goto L10;
L200:
	at1=0.;
	at2=-g;
	rt1d=vt1;
	rt2d=vt2;
	vt1d=at1;
	vt2d=at2;
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
	
