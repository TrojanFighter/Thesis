#include <stdio.h>
#include <stdlib.h>
#include <math.h>
void lambert(double xic,double yic,double tfdes,double xf,double yf,
	double *vrx,double *vry,double xlongm,double xlongt);
main()
{
	double xlongmdeg,xlongtdeg,altnmt,altnmm,tf,pi,degrad,a,gm;
	double altt,altm,xlongm,xlongt,xm,ym,xt,yt,vrxm,vrym;
	FILE *fptr;
	fptr=fopen("DATFIL.TXT","w");
	xlongmdeg=45.;
	xlongtdeg=90.;
	altnmt=0.;
	altnmm=0.;
	tf=1000.;
	pi=3.14159;
	degrad=360./(2.*pi);
	a=2.0926e7;
	gm=1.4077e16;
	altt=altnmt*6076.;
	altm=altnmm*6076.;
	xlongm=xlongmdeg/degrad;
	xlongt=xlongtdeg/degrad;
	xm=(a+altm)*cos(xlongm);
	ym=(a+altm)*sin(xlongm);
	xt=(a+altt)*cos(xlongt);
	yt=(a+altt)*sin(xlongt);
	lambert(xm,ym,tf,xt,yt,&vrxm,&vrym,xlongm,xlongt);
	printf("%10.4f %10.4f %10.4f\n",tf,vrxm,vrym);
	fprintf(fptr,"%10.4f %10.4f %10.4f\n",tf,vrxm,vrym);
	fclose (fptr);
	return 0;
}
	
	
void lambert(double xic,double yic,double tfdes,double xf,double yf,
	double *vrx,double *vry,double xlongm,double xlongt)
{
	int icount;
	double a,gm,ric,rf,cphi,phi,r0,pi,degrad,gamdeg,gam,top,temp,bot;
	double v,xlam,top1,bot1p,bot1,top2,bot2;
	double top3,bot3,tf,gamdegnew,gamdegfin;
	a=2.0926e7;
	gm=1.4077e16;
	ric=sqrt(xic*xic+yic*yic);
	rf=sqrt(xf*xf+yf*yf);
	cphi=(xic*xf+yic*yf)/(ric*rf);
	phi=atan2(sqrt(1.-cphi*cphi),cphi);
	r0=ric;
	pi=3.14159;
	degrad=360./(2.*pi);
	icount=0;
	for (gamdeg=-90.; gamdeg<=90.; gamdeg=gamdeg+.1){
		gam=gamdeg/degrad;
		top=gm*(1.-cos(phi));
		temp=r0*cos(gam)/rf-cos(phi+gam);
		bot=r0*cos(gam)*temp;
		if(top<0.||bot<0.)goto L10;
		v=sqrt(top/bot);
		if (xlongt>xlongm){
	  		*vrx=v*cos(pi/2.-gam+xlongm);
	  		*vry=v*sin(pi/2.-gam+xlongm);
	  	}
		else{
	  		*vrx=v*cos(-pi/2.+gam+xlongm);
	  		*vry=v*sin(-pi/2.+gam+xlongm);
		}
		xlam=r0*v*v/gm;
		top1=tan(gam)*(1-cos(phi))+(1-xlam)*sin(phi);
		bot1p=(1-cos(phi))/(xlam*cos(gam)*cos(gam));
		bot1=(2-xlam)*(bot1p+cos(gam+phi)/cos(gam));
		top2=2*cos(gam);
		if((2./xlam-1.)<0.)goto L10;
		bot2=xlam*pow((2./xlam-1.),1.5);
		top3=sqrt(2./xlam-1);
		bot3=cos(gam)/tan(phi/2)-sin(gam);
		temp=(top2/bot2)*atan2(top3,bot3);
		tf=r0*(top1/bot1+temp)/(v*cos(gam));
		if(tf>tfdes)goto L11;
L10:
	;
	}
L11:
 	gamdegnew=gamdeg-.15;
	gamdegfin=gamdeg+1.;
	for (gamdeg=gamdegnew; gamdeg<=gamdegfin; gamdeg=gamdeg+.0001){
		gam=gamdeg/degrad;
		top=gm*(1.-cos(phi));
		temp=r0*cos(gam)/rf-cos(phi+gam);
		bot=r0*cos(gam)*temp;
		if(top<0.||bot<0.)goto L20;
		v=sqrt(top/bot);
		if (xlongt>xlongm){
	  		*vrx=v*cos(pi/2.-gam+xlongm);
	  		*vry=v*sin(pi/2.-gam+xlongm);
	  	}
		else{
	  		*vrx=v*cos(-pi/2.+gam+xlongm);
	  		*vry=v*sin(-pi/2.+gam+xlongm);
		}
		xlam=r0*v*v/gm;
		top1=tan(gam)*(1-cos(phi))+(1-xlam)*sin(phi);
		bot1p=(1-cos(phi))/(xlam*cos(gam)*cos(gam));
		bot1=(2-xlam)*(bot1p+cos(gam+phi)/cos(gam));
		top2=2*cos(gam);
		if((2./xlam-1.)<0.)goto L20;
		bot2=xlam*pow((2./xlam-1.),1.5);
		top3=sqrt(2/xlam-1);
		bot3=cos(gam)/tan(phi/2)-sin(gam);
		temp=(top2/bot2)*atan2(top3,bot3);
		tf=r0*(top1/bot1+temp)/(v*cos(gam));
		if(tf>tfdes)goto L21;
 L20:
 	;
 	}
 L21:
 	;
 }
