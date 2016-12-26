#include <stdio.h>
#include <stdlib.h>
#include <math.h>
void lambert(double xm,double ym,double tf,double xtf,double ytf,
	double *vrxm,double *vrym,double xlongm,double xlongt);
void distance(double xt,double yt,double xfirstt,double yfirstt,double *distnmt);
void predict (double t,double tf,double xt,double yt,double x1t,double y1t,
		double *xtf,double *ytf,double wp1,double wtot,double tb1,
		double trst1,double tb2,double wp2,double wtot2,double trst2,
		double wpay,double xm,double ym,double x1m,double y1m,
		double *zem1,double *zem2,double tgo);
main()
{
	int step,boost,left;
	double h,a,gm,xnp,axmguid,aymguid,prederr,xisp1,xisp2,xmf1,xmf2;
	double wpay,delv,delv1,delv2,amax1,amax2,xkickdeg,top2,bot2,wp2;
	double ws2,wtot2,trst2,tb2,top1,bot1,wp1,ws1,wtot,trst1,tb1;
	double altnmt,altnmm,altt,altm,pi,degrad,s,scount,xlongmdeg;
	double xlongtdeg,xlongm,xlongt,xm,ym,xt,yt,xfirstt,yfirstt,x1t,y1t;
	double x1m,y1m,rtm1,rtm2,vtm1,vtm2,vc;
	double tgo,xoldt,yoldt,x1oldt,y1oldt,xoldm,yoldm,x1oldm,y1oldm;
	double xdt,ydt,x1dt,y1dt,xdm,ydm,x1dm,y1dm,delvd;
	double zemplos,zem1,zem2,xtf,ytf,xnc,xnmt,ynmt,xnmm,ynmm;
	double distnmt,distnmm,at,vel,axt,ayt,atplos,tembott,atplosg,xncg;
	double xlam,xlamd,tembotm;
	double t,tf,vrxm,vrym,rtm,delvold,wgt,trst;
	FILE *fptr;
	fptr=fopen("DATFIL.TXT","w");
	left=0;
 	h=.01;
	a=2.0926e7;
	gm=1.4077e16;
	xnp=3.;
	axmguid=0.;
	aymguid=0.;
	prederr=0.;
	xisp1=250.;
	xisp2=250.;
	xmf1=.85;;
	xmf2=.85;;
	wpay=100.;
	delv=20000.;
	delv1=.3333*delv;
	delv2=.6667*delv;
	amax1=20.;
	amax2=20.;
	xkickdeg=80.;
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
	altnmt=0.;
	altnmm=0.;
	altt=altnmt*6076.;
	altm=altnmm*6076.;
	pi=3.14159;
	degrad=360./(2.*pi);
	s=0.;
	scount=0.;
	xlongmdeg=85.;
	xlongtdeg=90.;
	xlongm=xlongmdeg/degrad;
	xlongt=xlongtdeg/degrad;
	xm=(a+altm)*cos(xlongm);
	ym=(a+altm)*sin(xlongm);
	xt=(a+altt)*cos(xlongt);
	yt=(a+altt)*sin(xlongt);
	xfirstt=xt;
	yfirstt=yt;
	if (left==1){
	  x1t=cos(pi/2.-xkickdeg/degrad+xlongt);
	  y1t=sin(pi/2.-xkickdeg/degrad+xlongt);
	}
	else{
	  x1t=cos(-pi/2.+xkickdeg/degrad+xlongt);
	  y1t=sin(-pi/2.+xkickdeg/degrad+xlongt);
	}
	t=0.;
	tf=50.;
	tgo=tf-t;
	predict (t,tf,xt,yt,x1t,y1t,&xtf,&ytf,wp1,wtot,tb1,
		trst1,tb2,wp2,wtot2,trst2,wpay,xm,ym,x1m,y1m,
		&zem1,&zem2,tgo);
	ytf=ytf+prederr;
	lambert(xm,ym,tf,xtf,ytf,&vrxm,&vrym,xlongm,xlongt);
	x1m=vrxm;
	y1m=vrym;
	rtm1=xt-xm;
	rtm2=yt-ym;
	rtm=sqrt(rtm1*rtm1+rtm2*rtm2);
	vtm1=x1t-x1m;
	vtm2=y1t-y1m;
	vc=-(rtm1*vtm1+rtm2*vtm2)/rtm;
	delv=0.;
L10:	if(vc<0.)goto L999;
 	tgo=rtm/vc;
	if(tgo>.1)
		h=.01;
	else
		h=.0001;
	xoldt=xt;
	yoldt=yt;
	x1oldt=x1t;
	y1oldt=y1t;
	xoldm=xm;
	yoldm=ym;
	x1oldm=x1m;
	y1oldm=y1m;
	delvold=delv;
	step=1;
	goto L200;
L66:	step=2;
	xt=xt+h*xdt;
	yt=yt+h*ydt;
	x1t=x1t+h*x1dt;
	y1t=y1t+h*y1dt;
	xm=xm+h*xdm;
	ym=ym+h*ydm;
	x1m=x1m+h*x1dm;
	y1m=y1m+h*y1dm;
	delv=delv+h*delvd;
	t=t+h;
	goto L200;
L55:	xt=(xoldt+xt)/2+.5*h*xdt;
	yt=(yoldt+yt)/2+.5*h*ydt;
	x1t=(x1oldt+x1t)/2+.5*h*x1dt;
	y1t=(y1oldt+y1t)/2+.5*h*y1dt;
	xm=(xoldm+xm)/2+.5*h*xdm;
	ym=(yoldm+ym)/2+.5*h*ydm;
	x1m=(x1oldm+x1m)/2+.5*h*x1dm;
	y1m=(y1oldm+y1m)/2+.5*h*y1dm;
	delv=(delvold+delv)/2.+.5*h*delvd;
	altt=sqrt(xt*xt+yt*yt)-a;
	altm=sqrt(xm*xm+ym*ym)-a;
	s=s+h;
	scount=scount+h;
	if(scount<.99999)goto L10;
 	scount=0.;
 	predict (t,tf,xt,yt,x1t,y1t,&xtf,&ytf,wp1,wtot,tb1,
		trst1,tb2,wp2,wtot2,trst2,wpay,xm,ym,x1m,y1m,
		&zem1,&zem2,tgo);
	zemplos=-zem1*sin(xlam)+zem2*cos(xlam);
	xnc=xnp*zemplos/(tgo*tgo);
	if(xnc>966.)xnc=966;
	if(xnc<-966.)xnc=-966.;
	xnmt=xt/6076.;
	ynmt=yt/6076.;
	xnmm=xm/6076.;
	ynmm=ym/6076.;
	altnmt=altt/6076.;
	distance(xt,yt,xfirstt,yfirstt,&distnmt);
	altnmm=altm/6076.;
	distance(xm,ym,xfirstt,yfirstt,&distnmm);
	xncg=xnc/32.2;
	atplosg=atplos/32.2;
	printf("%10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n",t,distnmt
			,altnmt,distnmm,altnmm,xncg,delv,atplosg);
	fprintf(fptr,"%10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n",t,distnmt
			,altnmt,distnmm,altnmm,xncg,delv,atplosg);
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
	vel=sqrt(x1t*x1t+y1t*y1t);
	axt=at*x1t/vel;
	ayt=at*y1t/vel;
	atplos=ayt*cos(xlam)-axt*sin(xlam);
	tembott=pow((xt*xt+yt*yt),1.5);
	x1dt=-gm*xt/tembott+axt;
	y1dt=-gm*yt/tembott+ayt;
	xdt=x1t;
	ydt=y1t;
	rtm1=xt-xm;
	rtm2=yt-ym;
	rtm=sqrt(rtm1*rtm1+rtm2*rtm2);
	vtm1=x1t-x1m;
	vtm2=y1t-y1m;
	vc=-(rtm1*vtm1+rtm2*vtm2)/rtm;
	tgo=rtm/vc;
	xlam=atan2(rtm2,rtm1);
	xlamd=(rtm1*vtm2-rtm2*vtm1)/(rtm*rtm);
	axmguid=-xnc*sin(xlam);
	aymguid=xnc*cos(xlam);
	delvd=fabs(xnc);
	tembotm=pow((xm*xm+ym*ym),1.5);
	x1dm=-gm*xm/tembotm+axmguid;
	y1dm=-gm*ym/tembotm+aymguid;
	xdm=x1m;
	ydm=y1m;
	if(step<=1)
		goto L66;
	else
		goto L55;
L999:
 	xnmt=xt/6076.;
	ynmt=yt/6076.;
	xnmm=xm/6076.;
	ynmm=ym/6076.;
	altnmt=altt/6076.;
	distance(xt,yt,xfirstt,yfirstt,&distnmt);
	altnmm=altm/6076.;
	distance(xm,ym,xfirstt,yfirstt,&distnmm);
	xncg=xnc/32.2;
	atplosg=atplos/32.2;
	printf("%10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n",t,distnmt
			,altnmt,distnmm,altnmm,xncg,delv,atplosg);
	fprintf(fptr,"%10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n",t,distnmt
			,altnmt,distnmm,altnmm,xncg,delv,atplosg);
	printf("%10.4f %10.4f %10.4f\n",t,rtm,delv);
 	fclose (fptr);
 	return 0;
}
	
void distance(double xt,double yt,double xfirstt,double yfirstt,double *distnmt)
{
	double r,a,cbeta,beta;
	r=sqrt(xt*xt+yt*yt);
	a=2.0926E7;
	cbeta=(xt*xfirstt+yt*yfirstt)/(r*a);
	if(cbeta<1.){
		beta=atan2(sqrt(1.-cbeta*cbeta),cbeta);
		*distnmt=a*beta/6076.;
	}
	else
		*distnmt=(xfirstt-xt)/6076.;
}

void predict (double t,double tf,double xt,double yt,double x1t,double y1t,
		double *xtf,double *ytf,double wp1,double wtot,double tb1,
		double trst1,double tb2,double wp2,double wtot2,double trst2,
		double wpay,double xm,double ym,double x1m,double y1m,
		double *zem1,double *zem2,double tgo)
{
	int step;
	double h,a,gm,x,y,x1,y1xold,yold,x1old,y1old,xoldm,yoldm;
	double x1oldm,y1oldm,xd,yd,x1d,y1d,xdm,ydm,x1dm,y1dm;
	double wgt,trst,at,vel,axt,ayt,tembott,tembotm,y1,xold;
	if(tgo>1)
		h=.01;
	else
		h=tgo;;
	a=2.0926e7;
	gm=1.4077e16;
	x=xt;
	y=yt;
	x1=x1t;
	y1=y1t;
L10:	if(t>(tf-.00001))goto L999;
	xold=x;
	yold=y;
	x1old=x1;
	y1old=y1;
	xoldm=xm;
	yoldm=ym;
	x1oldm=x1m;
	y1oldm=y1m;
	step=1;
	goto L200;
L66:	step=2;
	x=x+h*xd;
	y=y+h*yd;
	x1=x1+h*x1d;
	y1=y1+h*y1d;
	xm=xm+h*xdm;
	ym=ym+h*ydm;
	x1m=x1m+h*x1dm;
	y1m=y1m+h*y1dm;
	t=t+h;
	goto L200;
L55:	x=(xold+x)/2+.5*h*xd;
	y=(yold+y)/2+.5*h*yd;
	x1=(x1old+x1)/2+.5*h*x1d;
	y1=(y1old+y1)/2+.5*h*y1d;
	xm=(xoldm+xm)/2+.5*h*xdm;
	ym=(yoldm+ym)/2+.5*h*ydm;
	x1m=(x1oldm+x1m)/2+.5*h*x1dm;
	y1m=(y1oldm+y1m)/2+.5*h*y1dm;
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
	tembott=pow((x*x+y*y),1.5);
	x1d=-gm*x/tembott+axt;
	y1d=-gm*y/tembott+ayt;
	xd=x1;
	yd=y1;
	tembotm=pow((xm*xm+ym*ym),1.5);
	x1dm=-gm*xm/tembotm;
	y1dm=-gm*ym/tembotm;
	xdm=x1m;
	ydm=y1m;
	if(step<=1)
		goto L66;
	else
		goto L55;
L999:
 	*xtf=x;
	*ytf=y;
	*zem1=x-xm;
	*zem2=y-ym;
}
	
void lambert(double xm,double ym,double tfdes,double xt,double yt,
		double*vrx,double*vry,double xlongm,double xlongt)
{
	double a,gm,ric,rf,cphi,phi,sphi,r0,pi,degrad,gmin,gmax,gam;
	int icount;
	double top,temp,bot,v,xlam,top1,bot1p,bot1,top2,bot2;
	double top3,bot3,tf,xnext,told,gold;
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
	;
}
