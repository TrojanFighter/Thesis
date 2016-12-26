#include <stdio.h>
#include <stdlib.h>
#include <math.h>
void lambert(double x,double y,double tgolam,double xf,double yf,
	double *vrx,double *vry,double xlongm,double xlongt);
void distance(double x,double y,double xfirst,double yfirst,double *distnm);
main()
{
	int step,left,qboost,qzero;
	double xisp1,xisp2,xmf1,xmf2,wpay,delv,delv1,delv2,amax1;
	double amax2,top2,bot2,wp2,ws2,wtot2,trst2,tb2,top1,bot1,wp1,ws1;
	double wtot,trst1,tb1,delvk,h,t,s,a,gm,altnm,alt,angdeg,ang,xlongm;
	double x,y,x1,y1,axt,ayt,xlongtdeg,xlongtxf,yf,xfirst,yfirst;
	double tf,dvcap,xold,yold,x1old,y1old,xd,yd,x1d,y1d;
	double tgolam,at,delx,dely,del,thet,degthet,phi,degphi,xlongt;
	double velk,xnm,ynm,wgt,trst,tembot,xf,distnm,distinitnm,vrx,vry;
	FILE *fptr;
	fptr=fopen("DATFIL.TXT","w");
	left=1;
	qboost=1;
	qzero=0;
	xisp1=300.;
	xisp2=300.;
	xmf1=.90;
	xmf2=.90;
	wpay=100.;
	delv=20000.;
	delv1=.3333*delv;
	delv2=.6667*delv;
	amax1=20.;
	amax2=20.;
	top2=wpay*(exp(delv2/(xisp2*32.2))-1.);
	bot2=1./xmf2-((1.-xmf2)/xmf2)*exp(delv2/(xisp2*32.2));
	wp2=top2/bot2;
	ws2=wp2*(1.-xmf2)/xmf2;
	wtot2=wp2+ws2+wpay;
	trst2=amax2*(wpay+ws2);
	tb2=xisp2*wp2/trst2;
	top1=wtot2*(exp(delv1/(xisp1*32.2))-1.);
	bot1=1./xmf1-((1.-xmf1)/xmf1)*exp(delv1/(xisp1*32.2));
	wp1=top1/bot1;
	ws1=wp1*(1.-xmf1)/xmf1;
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
	angdeg=30.;
	ang=angdeg/57.3;
	xlongm=ang;
	x=(a+alt)*cos(ang);
	y=(a+alt)*sin(ang);
	alt=sqrt(x*x+y*y)-a;
	x1=0.;
	y1=0.;
	axt=0.;
	ayt=0.;
	xlongtdeg=45.;
	xlongt=xlongtdeg/57.3;
	xf=a*cos(xlongt);
	yf=a*sin(xlongt);
	xfirst=xf;
	yfirst=yf;
	distance(x,y,xfirst,yfirst,&distnm);
	distinitnm=distnm;
	tf=500.;
	dvcap=delv;
L10:
	if(alt<0.&&t>10.)
		goto L999;
	xold=x;
	yold=y;
	x1old=x1;
	y1old=y1;
	step=1;
	goto L200;
L66:	
	step=2;
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
	tgolam=tf-t;
	dvcap=dvcap-h*at;
	if(qboost==1&&dvcap>50.){
		xlongm=atan2(y,x);
		lambert(x,y,tgolam,xf,yf,&vrx,&vry,xlongm,xlongt);
		delx=vrx-x1;
		dely=vry-y1;
		del=sqrt(delx*delx+dely*dely);
		if(qzero==0&&dvcap>del){
			thet=sqrt(6.*(1.-del/dvcap));
			degthet=57.3*thet;
		}
		else
			qzero=1;
		phi=atan2(dely,delx);
		degphi=57.3*phi;
		if(xlongt>xlongm){
			axt=at*cos(phi-thet);
			ayt=at*sin(phi-thet);
		}
		else{
			axt=at*cos(phi+thet);
			ayt=at*sin(phi+thet);
		}
		distance(x,y,xfirst,yfirst,&distnm);
		distnm=distinitnm-distnm;
		altnm=(sqrt(x*x+y*y)-a)/6076.;
	}
	else if(qboost==1){
		lambert(x,y,tgolam,xf,yf,&vrx,&vry,xlongm,xlongt);
		trst=0.;
		qboost=0;
		axt=0.;
		ayt=0.;
		x1=vrx;
		y1=vry;
		x1old=x1;
		y1old=y1;
		distance(x,y,xfirst,yfirst,&distnm);
		distnm=distinitnm-distnm;
		altnm=(sqrt(x*x+y*y)-a)/6076.;
	}
	else{
		qboost=0;
		axt=0.;
		ayt=0.;
	}
	if(s<9.99999)
		goto L10;
 	s=0.;
	distance(x,y,xfirst,yfirst,&distnm);
	distnm=distinitnm-distnm;
	altnm=(sqrt(x*x+y*y)-a)/6076.;
	velk=sqrt(x1*x1+y1*y1)/1000.;
	xnm=x/6076.;
	ynm=y/6076.;
	printf("%10.4f %10.4f %10.4f\n",t,distnm,altnm);
	fprintf(fptr,"%10.4f %10.4f %10.4f\n",t,distnm,altnm);
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
	xd=x1;
	yd=y1;
	tembot=pow((x*x+y*y),1.5);
	x1d=-gm*x/tembot+axt;
	y1d=-gm*y/tembot+ayt;
	if(step<=1)
		goto L66;
	else
		goto L55;
L999:
 	distance(x,y,xfirst,yfirst,&distnm);
	distnm=distinitnm-distnm;
	altnm=(sqrt(x*x+y*y)-a)/6076.;
	velk=sqrt(x1*x1+y1*y1)/1000.;
	xnm=x/6076.;
	ynm=y/6076.;
	printf("%10.4f %10.4f %10.4f\n",t,distnm,altnm);
	fprintf(fptr,"%10.4f %10.4f %10.4f\n",t,distnm,altnm);
	fclose (fptr);
	return 0;
}
	
	void distance(double x,double y,double xfirst,double yfirst,double *distnm)
{
	double r,a,cbeta,beta;
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
