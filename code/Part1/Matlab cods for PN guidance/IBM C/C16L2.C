#include <stdio.h>
#include <stdlib.h>
#include <math.h>
void lambert(double xm,double ym,double tf,double xtf,double ytf,
	double *vrxm,double *vrym,double xlongm,double xlongt);
void distance(double xt,double yt,double xfirstt,double yfirstt,double *distnmt);
void gauss(double siglam,double *xlamnoise);
void predict (double tf,double xt,double yt,double x1t,double y1t,
		double *xtf,double *ytf,double wp1,double wtot,double tb1,
		double trst1,double tb2,double wp2,double wtot2,
		double trst2,double wpay);
main()
{
	int step;
	double h,a,gm,xnp,ts,xlongtdeg,xlongmdeg,degrad,siglam,prederr;
	double xisp1,xisp2,xmf1,xmf2,wpay,delv,delv1,delv2,amax1,amax2;
	double xkickdeg,top2,bot2,wp2,ws2,wtot2,trst2,tb2,top1,bot1;
	double wp1,ws1,wtot,trst1,tb1,altnmt,altnmm,altt,altm,s,scount;
	double xlongm,xlongt,xm,ym,xt,yt,xfirstt,yfirstt,x1t,y1t;
	double k1,k2,k3,xth,xtdh,xtddh,yth,ytdh,ytddh,t,tf,phin;
	double xtf,ytf,x1m,y1m,vrxm,vrym,rtm1,rtm2,rtm,sigpos;
	double p11,p12,p13,p22,p23,p33,vtm1,vtm2,vc,tgo;
	double xoldt,yoldt,x1oldt,y1oldt,xoldm,yoldm,x1oldm,y1oldm,delvold;
	double xdt,ydt,x1dt,y1dt,xdm,ydm,x1dm,y1dm,delvd;
	double ts2,ts3,ts4,ts5,sign2,m11,m12,m13,m22,m23,m33,bot,fact;
	double xlamnoise,ytmeas,xtmeas,resx,resy,xmnt,ymnt,xnmm,ynmm;
	double distnmt,distnmm,trst,wgt,at,vel,axt,ayt,tembott;
	double atplos,xlam,xlamd,xnc,am1,am2,tembotm,xnmt,ynmt;
	double k1p,k2p,k3p,sigx2,sigx,sigy2,sigy;
	double p11p,p12p,p13p,p22p,p23p,p33p;
	double m11p,m12p,m13p,m22p,m23p,m33p,botp,factp;
	double x1dtg,xtddhg,y1dtg,ytddhg;
	FILE *fptr;
	fptr=fopen("DATFIL.TXT","w");
	a=2.0926e7;
	gm=1.4077e16;
	xnp=3.;
	ts=1.;
	xlongtdeg=90.;
	xlongmdeg=85.;
	degrad=57.3;
	siglam=.001;
	prederr=0.;
	xisp1=250.;
	xisp2=250.;
	xmf1=.85;
	xmf2=.85;
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
	s=0.;
	scount=0.;
	xlongm=xlongmdeg/degrad;
	xlongt=xlongtdeg/degrad;
	xm=a*cos(xlongm);
	ym=a*sin(xlongm);
	xt=(a+altt)*cos(xlongt);
	yt=(a+altt)*sin(xlongt);
	xfirstt=xt;
	yfirstt=yt;
	x1t=cos(xkickdeg/57.3);
	y1t=sin(xkickdeg/57.3);
	k1=0.;
	k2=0.;
	k3=0.;
	k1p=0.;
	k2p=0.;
	k3p=0.;
	xth=xt;
	xtdh=x1t;
	xtddh=0.;
	yth=yt;
	ytdh=y1t;
	ytddh=0.;
	t=0.;
	tf=50.;
	phin=100*100/tf;
	predict (tf,xt,yt,x1t,y1t,&xtf,&ytf,wp1,wtot,tb1,
		trst1,tb2,wp2,wtot2,trst2,wpay);
	ytf=ytf+prederr;
	lambert(xm,ym,tf,xtf,ytf,&vrxm,&vrym,xlongm,xlongt);
	x1m=vrxm;
	y1m=vrym;
	rtm1=xt-xm;
	rtm2=yt-ym;
	rtm=sqrt(rtm1*rtm1+rtm2*rtm2);
	xlam=atan2(rtm2,rtm1);
	sigx2=(siglam*rtm*sin(xlam))*(siglam*rtm*sin(xlam));
	sigx=sqrt(sigx2);
	p11=sigx2;
	p12=0.;
	p13=0.;
	p22=0.;
	p23=0.;
	p33=100*100;
	sigy2=(siglam*rtm*cos(xlam))*(siglam*rtm*cos(xlam));
	sigy=sqrt(sigy2);
	p11p=sigy2;
	p12p=0.;
	p13p=0.;
	p22p=0.;
	p23p=0.;
	p33p=100*100;
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
	if(scount<(ts-.00001))goto L10;
 	scount=0.;
	ts2=ts*ts;
	ts3=ts2*ts;
	ts4=ts3*ts;
	ts5=ts4*ts;
	sigx2=(siglam*rtm*sin(xlam))*(siglam*rtm*sin(xlam));
	sigy2=(siglam*rtm*cos(xlam))*(siglam*rtm*cos(xlam));
	sigx=sqrt(sigx2);
	sigy=sqrt(sigy2);
	m11=p11+ts*p12+.5*ts2*p13+ts*(p12+ts*p22+.5*ts2*p23);
	m11=m11+.5*ts2*(p13+ts*p23+.5*ts2*p33)+ts5*phin/20.;
	m12=p12+ts*p22+.5*ts2*p23+ts*(p13+ts*p23+.5*ts2*p33)+ts4*phin/8.;
	m13=p13+ts*p23+.5*ts2*p33+phin*ts3/6.;
	m22=p22+ts*p23+ts*(p23+ts*p33)+phin*ts3/3.;
	m23=p23+ts*p33+.5*ts2*phin;
	m33=p33+phin*ts;
	bot=m11+sigx2;
	k1=m11/bot;
	k2=m12/bot;
	k3=m13/bot;
	fact=1.-k1;
	p11=fact*m11;
	p12=fact*m12;
	p13=fact*m13;
	p22=-k2*m12+m22;
	p23=-k2*m13+m23;
	p33=-k3*m13+m33;
	m11p=p11p+ts*p12p+.5*ts2*p13p+ts*(p12p+ts*p22p+.5*ts2*p23p);
	m11p=m11p+.5*ts2*(p13p+ts*p23p+.5*ts2*p33p)+ts5*phin/20.;
	m12p=p12p+ts*p22p+.5*ts2*p23p+ts*(p13p+ts*p23p+.5*ts2*p33p)+
     		ts4*phin/8.;
	m13p=p13p+ts*p23p+.5*ts2*p33p+phin*ts3/6.;
	m22p=p22p+ts*p23p+ts*(p23p+ts*p33p)+phin*ts3/3.;
	m23p=p23p+ts*p33p+.5*ts2*phin;
	m33p=p33p+phin*ts;
	botp=m11p+sigy2;
	k1p=m11p/botp;
	k2p=m12p/botp;
	k3p=m13p/botp;
	factp=1.-k1p;
	p11p=factp*m11p;
	p12p=factp*m12p;
	p13p=factp*m13p;
	p22p=-k2p*m12p+m22p;
	p23p=-k2p*m13p+m23p;
	p33p=-k3p*m13p+m33p;
	gauss(siglam,&xlamnoise);
	ytmeas=ym+rtm*sin(xlam+xlamnoise);
	xtmeas=xm+rtm*cos(xlam+xlamnoise);
	resx=xtmeas-xth-ts*xtdh-.5*ts2*xtddh;
	xth=k1*resx+xth+ts*xtdh+.5*ts2*xtddh;
	xtdh=k2*resx+xtdh+ts*xtddh;
	xtddh=k3*resx+xtddh;
	resy=ytmeas-yth-ts*ytdh-.5*ts2*ytddh;
	yth=k1p*resy+yth+ts*ytdh+.5*ts2*ytddh;
	ytdh=k2p*resy+ytdh+ts*ytddh;
	ytddh=k3p*resy+ytddh;
	xnmt=xt/6076.;
	ynmt=yt/6076.;
	xnmm=xm/6076.;
	ynmm=ym/6076.;
	altnmt=altt/6076.;
	distance(xt,yt,xfirstt,yfirstt,&distnmt);
	altnmm=altm/6076.;
	distance(xm,ym,xfirstt,yfirstt,&distnmm);
	x1dtg=x1dt/32.2;
	xtddhg=xtddh/32.2;
	y1dtg=y1dt/32.2;
	ytddhg=ytddh/32.2;
	printf("%10.3f %10.3f %10.3f %10.3f %10.3f \n",t,x1dtg,xtddhg,y1dtg,ytddhg);
	fprintf(fptr,"%10.3f %10.3f %10.3f %10.3f %10.3f \n",t,x1dtg,xtddhg,y1dtg,ytddhg);
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
	tembott=pow((xt*xt+yt*yt),1.5);
	x1dt=-gm*xt/tembott+axt;
	y1dt=-gm*yt/tembott+ayt;
	atplos=y1dt*cos(xlam)-x1dt*sin(xlam);
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
	xnc=xnp*vc*xlamd+.5*xnp*atplos;
	delvd=fabs(xnc);
	am1=-xnc*sin(xlam);
	am2=xnc*cos(xlam);
	tembotm=pow((xm*xm+ym*ym),1.5);
	x1dm=-gm*xm/tembotm+am1;
	y1dm=-gm*ym/tembotm+am2;
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

void predict (double tf,double xt,double yt,double x1t,double y1t,
		double *xtf,double *ytf,double wp1,double wtot,double tb1,
		double trst1,double tb2,double wp2,double wtot2,
		double trst2,double wpay)
{
	int step;
	double h,a,gm,t,x,y,x1,y1,xold,yold,x1old,y1old,xd,yd,x1d,y1d;
	double wgt,trst,at,vel,axt,ayt,tembott;
	h=.01;
	a=2.0926e7;
	gm=1.4077e16;
	t=0.;
	x=xt;
	y=yt;
	x1=x1t;
	y1=y1t;
L10:	if(t>(tf-.00001))goto L999;
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
	if(step<=1)
		goto L66;
	else
		goto L55;
L999:
 	*xtf=x;
	*ytf=y;
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
	
void gauss(double signoise,double *xlamnoise)
{
	double sum,x,y,temp;
	int j;
	sum=0.;
	for (j=1; j<=6; j=j+1){
		x=rand();
		y=(float)x/(float)RAND_MAX;
		sum=sum+y;
		}
	temp=sum-3.;
	*xlamnoise=1.414*temp*signoise;
}
