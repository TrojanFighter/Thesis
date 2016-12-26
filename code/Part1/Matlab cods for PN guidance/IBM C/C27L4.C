#include <stdio.h>
#include <stdlib.h>
#include <math.h>
void predict(double t,double xm,double ym,double xmd,double ymd,double phi,
		double tf,double phiad,double phibd,double *x1,double *y1,
		double xnc,double h,double phidmax,double accerr);
main()
{
	int step;
	double gain,xt,yt,phiaddeg,phibddeg,phidmaxdeg,tf,accerr,h,ts,xnc;
	double phidmax,phi,xm,ym,xmd,ymd,phiad,phibd,s,t,xmold,ymold,xmdold;
	double ymdold,phiold,delx,dely,phibd2,dxdpb,dydpb,phiad3,x1,y1,x3,y3,x2,y2;
	double dxdpa,dydpa,delphiad,delphibd,slope,bint,phid,xmdd,ymdd,rtm;
	FILE *fptr1;
	fptr1=fopen("DATFIL.TXT","w");
	gain=1.;
	xt=-4000.;
	yt=5000.;
	phiaddeg=10.;
	phibddeg=5.;
	phidmaxdeg=15.;
	tf=100.;
	accerr=0.;
	h=.01;
	ts=.1;
	xnc=12.;
	phidmax=phidmaxdeg/57.3;
	phi=0.;
	xm=0.;
	ym=0.;
	xmd=0.;
	ymd=0.;
	phiad=phiaddeg/57.3;
	phibd=phibddeg/57.3;
	s=0.;
	t=0.;
L10:
	if(t>(tf-.00001))goto L999;
	xmold=xm;
	ymold=ym;
	xmdold=xmd;
	ymdold=ymd;
	phiold=phi;
	step=1;
	goto L200;
L66:
	step=2;
	xm=xm+h*xmd;
	ym=ym+h*ymd;
	xmd=xmd+h*xmdd;
	ymd=ymd+h*ymdd;
	phi=phi+h*phid;
	t=t+h;
	goto L200;
L55:
	xm=(xmold+xm)/2+.5*h*xmd;
	ym=(ymold+ym)/2+.5*h*ymd;
	xmd=(xmdold+xmd)/2+.5*h*xmdd;
	ymd=(ymdold+ymd)/2+.5*h*ymdd;
	phi=(phiold+phi)/2.+.5*h*phid;
	s=s+h;
	if(s<(ts-.0001))goto L10;
 	s=0.;
	predict(t,xm,ym,xmd,ymd,phi,tf,phiad,phibd,&x1,&y1,xnc,h,phidmax,accerr);
	delx=xt-x1;
	dely=yt-y1;
	phibd2=phibd+.001;
	predict(t,xm,ym,xmd,ymd,phi,tf,phiad,phibd2,&x2,&y2,xnc,h,phidmax,accerr);
	dxdpb=(x2-x1)/(phibd2-phibd);
	dydpb=(y2-y1)/(phibd2-phibd);
	phiad3=phiad+.001;
	predict(t,xm,ym,xmd,ymd,phi,tf,phiad3,phibd,&x3,&y3,xnc,h,phidmax,accerr);
	dxdpa=(x3-x1)/(phiad3-phiad);
	dydpa=(y3-y1)/(phiad3-phiad);
	if((dxdpb*dydpa-dxdpa*dydpb)==0.){
		delphiad=0.;
		delphibd=0.;
	}
	else{
		delphiad=(dxdpb*dely-delx*dydpb)/(dxdpb*dydpa-dxdpa*dydpb);
		delphibd=(delx*dydpa-dxdpa*dely)/(dxdpb*dydpa-dxdpa*dydpb);
     	}
	phiad=phiad+gain*delphiad;
	phibd=phibd+gain*delphibd;
	printf("%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n"
		,t,xm,ym,xt,yt,phid*57.3);
	fprintf(fptr1,"%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n"
		,t,xm,ym,xt,yt,phid*57.3);
 	goto L10;
L200:
	slope=(phibd-phiad)/tf;
	bint=phibd-slope*tf;
	phid=slope*t+bint;
	if(phid>phidmax)phid=phidmax;
	if(phid<-phidmax)phid=-phidmax;
	xmdd=xnc*cos(phi);
	ymdd=xnc*sin(phi);
	if (step<2)
		goto L66;
	else
		goto L55;
L999:
 	rtm=sqrt((xt-xm)*(xt-xm)+(yt-ym)*(yt-ym));
 	printf("%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n"
		,t,xm,ym,xt,yt,phid*57.3);
	fprintf(fptr1,"%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n"
		,t,xm,ym,xt,yt,phid*57.3);
	printf("%6.3f\n",rtm);
	fclose (fptr1);
	return 0;
}
	
void predict(double tp,double xp,double yp,double xdp,double ydp,double thetp,
		double tf,double thadp,double thbdp,double *xf,double *yf,
     		double acc,double hp,double thdmax,double accerr)
{
	int step;
	double h,t,x,y,xd,yd,thet,thad,thbd,xold,yold,xdold,ydold,thetold;
	double slope,bint,thetd,xdd,ydd,s;
	h=hp*10.;
	t=tp;
	x=xp;
	y=yp;
	xd=xdp;
	yd=ydp;
	thet=thetp;
	thad=thadp;
	thbd=thbdp;
	s=0.;
L10:
	if(t>(tf-.00001))goto L999;
	xold=x;
	yold=y;
	xdold=xd;
	ydold=yd;
	thetold=thet;
	step=1;
	goto L200;
L66:
	step=2;
	x=x+h*xd;
	y=y+h*yd;
	xd=xd+h*xdd;
	yd=yd+h*ydd;
	thet=thet+h*thetd;
	t=t+h;
	goto L200;
L55:
	x=(xold+x)/2+.5*h*xd;
	y=(yold+y)/2+.5*h*yd;
	xd=(xdold+xd)/2+.5*h*xdd;
	yd=(ydold+yd)/2+.5*h*ydd;
	thet=(thetold+thet)/2.+.5*h*thetd;
	goto L10;
L200:
 	slope=(thbd-thad)/tf;
	bint=thbd-slope*tf;
	thetd=slope*t+bint;
	if(thetd>thdmax)thetd=thdmax;
	if(thetd<-thdmax)thetd=-thdmax;
	xdd=(acc+accerr)*cos(thet);
	ydd=(acc+accerr)*sin(thet);
	if (step<2)
		goto L66;
	else
		goto L55;
L999:
 	*xf=x;
	*yf=y;
 }

	
