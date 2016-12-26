#include <stdio.h>
#include <stdlib.h>
#include <math.h>
void lambert(double xm,double ym,double tfdes,double xt,double yt,
	double*vrx,double*vry,double xlongm,double xlongt);
main()
{
	double xlongmdeg,xlongtdeg,altnmt,altnmm,tfdes,pi,degrad,a,gm,altt;
	double altm,xlongm,xlongt,xm,ym,xt,yt,vrx,vry;
	xlongmdeg=45.;
	xlongtdeg=90.,
	altnmt=0.;
	altnmm=0.;
	tfdes=1000.;
	pi=3.14159;
	degrad=360./(2.*pi);
	a=2.0926E7;
	gm=1.4077E16;
	altt=altnmt*6076.;
	altm=altnmm*6076.;
	xlongm=xlongmdeg/degrad;
	xlongt=xlongtdeg/degrad;
	xm=(a+altm)*cos(xlongm);
	ym=(a+altm)*sin(xlongm);
	xt=(a+altt)*cos(xlongt);
	yt=(a+altt)*sin(xlongt);
	lambert(xm,ym,tfdes,xt,yt,&vrx,&vry,xlongm,xlongt);
	;
	return 0;
}
	void lambert(double xm,double ym,double tfdes,double xt,double yt,
		double*vrx,double*vry,double xlongm,double xlongt)
{
	FILE *fptr;
	double a,gm,ric,rf,cphi,phi,sphi,r0,pi,degrad,gmin,gmax,gam;
	int icount;
	double top,temp,bot,v,xlam,top1,bot1p,bot1,top2,bot2;
	double top3,bot3,tf,xnext,told,gold;
	fptr=fopen("DATFIL.TXT","w");
	a=2.0926E7;
	gm=1.4077E16;
	ric=sqrt(xm*xm+ym*ym);
	rf=sqrt(xt*xt+yt*yt);
	cphi=(xm*xt+ym*yt)/(ric*rf);
	phi=atan2(sqrt(1.-cphi*cphi),cphi);
	sphi=sin(phi);
	r0=ric;
	pi=3.14159;
	degrad=360./(2.*pi);
	gmin=atan2((sphi-sqrt(2.*r0*(1.-cphi)/rf)),(1.-cphi));
	gmax=atan2((sphi+sqrt(2.*r0*(1.-cphi)/rf)),(1.-cphi));
	gam=(gmin+gmax)/2.;
	
	for (icount=1; icount<=100; icount=icount+1){
		top=gm*(1.-cos(phi));
		temp=r0*cos(gam)/rf-cos(phi+gam);
		bot=r0*cos(gam)*temp;
		v=sqrt(top/bot);
		if(xlongt>xlongm){
	  		*vrx=v*cos(pi/2.-gam+xlongm);
	  		*vry=v*sin(pi/2.-gam+xlongm);
	  	}
		if(xlongt<=xlongm){
	  		*vrx=v*cos(-pi/2.+gam+xlongm);
	  		*vry=v*sin(-pi/2.+gam+xlongm);
	  	}
		xlam=r0*v*v/gm;
		top1=tan(gam)*(1.-cos(phi))+(1.-xlam)*sin(phi);
		bot1p=(1.-cos(phi))/(xlam*cos(gam)*cos(gam));
		bot1=(2.-xlam)*(bot1p+cos(gam+phi)/cos(gam));
		top2=2.*cos(gam);
		bot2=xlam*pow((2./xlam-1.),1.5);
		top3=sqrt(2./xlam-1.);
		bot3=cos(gam)/tan(phi/2.)-sin(gam);
		temp=(top2/bot2)*atan2(top3,bot3);
		tf=r0*(top1/bot1+temp)/(v*cos(gam));
		printf("%5d %10.4f %8.2f %8.2f %9.4f\n",icount,57.3*gam,*vrx,*vry,tf);
		fprintf(fptr,"%5d %10.4f %8.2f %8.2f %9.4f\n",icount,57.3*gam,*vrx,*vry,tf);
		if((fabs(tfdes-tf)<=.00000001*tfdes))
			goto L100;
		if(tf>tfdes)
			gmax=gam;
		else
			gmin=gam;
		
		if(icount==1)
			xnext=(gmax+gmin)/2.;
		else{
			xnext=gam+(gam-gold)*(tfdes-tf)/(tf-told);
			if(xnext>gmax||xnext<gmin)
				xnext=(gmax+gmin)/2.;
		}
		gold=gam;
		told=tf;
		gam=xnext;
	}
L100:
	fclose (fptr);
}
		
		
		
		
		
	

	
	
	
