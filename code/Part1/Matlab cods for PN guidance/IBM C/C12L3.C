#include <stdio.h>
#include <stdlib.h>
#include <math.h>
main()
{
	int step;
	float h,a,gm,gam,altnm,alt,xlam,v,angdeg,vrx,vry,s,scount;
	float x,y,xfirst,yfirst,x1,y1,t,tf,xold,yold,x1old,y1old;
	float xd,yd,x1d,y1d,xnm,ynm,tembot,ang;
	FILE *fptr;
	fptr=fopen("DATFIL.TXT","w");
 	h=.01;
	a=2.0926e7;
	gm=1.4077e16;
	gam=0.;
	altnm=1000.;
	alt=altnm*6076.;
	xlam=1.;
	v=sqrt(gm*xlam/(a+alt));
	angdeg=90.;
	ang=angdeg/57.3;
	vrx=v*cos(1.5708-gam/57.3+ang);
	vry=v*sin(1.5708-gam/57.3+ang);
	s=0.;
	scount=0.;
	x=(a+alt)*cos(ang);
	y=(a+alt)*sin(ang);
	xfirst=x;
	yfirst=y;
	x1=vrx;
	y1=vry;
	t=0.;
	tf=30000.;
L10:	if(t>tf)goto L999;
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
	s=s+h;
	scount=scount+h;
	if(scount<49.99999)goto L10;
 	scount=0.;
	xnm=x/6076.;
	ynm=y/6076.;
	printf("%8.2f %8.2f %8.2f\n",t,xnm,ynm);
	fprintf(fptr,"%8.2f %8.2f %8.2f\n",t,xnm,ynm);
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
	fclose (fptr);
	return 0;
}
