#include <stdio.h>
#include <stdlib.h>
#include <math.h>
void gauss(double signoise,double *xnoise);
void project(double ts,double yh,double ydh,double ytddh,double ytdddh,
	double *yb,double *ydb,double *ytddb,double *ytdddb,double hp,double xnl,
	double wh);
main()
{
	double phi[6][6],p[6][6],m[6][6],phip[6][6],phipphit[6][6],gain[6][2];
	double q[6][6],hmat[2][6],hm[2][6],mht[6][2];
	double phit[6][6],rmat[2][2],hmhtr[2][2],hmhtrinv[2][2];
	double hmht[2][2],ht[6][2],kh[6][6],idn[6][6],ikh[6][6];
	int order,step,apn,qperfect,i,j,k;
	double phis2,xnt,w,phasedeg,sigrin,siggl,srn,ra,whic,ts,tf,phis1;
	double vc,xnp,xnclim,tau,hedeg,vm,phase,x,tgo,t,s,y,yd,ytdd,ytddd;
	double xnc,xnl,h,hp,ts2,ts3,ts4,ts5,ts6,ts7,wh,yh,ydh,ytddh,ytdddh;
	double xnoise;
	double xlam,res,xs,top,bot1,bot2,xnpp,c1,c2,c3,c4,xp,c5,errytdd;
 	double errytddg,errytddd,errytdddg,sp55,sp55p,errw,sp44,sp44p;
 	double sp33,sp33p,sp33g,sp33pg,sp44g,sp44pg,ytddg,ytddhg,ytdddg,ytdddhg,xncg;
	double xnld,rtm,signoise,ynoise,yold,ydold,xnlold,ydd;
	double xlams,yb,ydb,ytddb,ytdddb;
	FILE *fptr1;
	FILE *fptr2;
	fptr1=fopen("DATFIL.TXT","w");
	fptr2=fopen("COVFIL.TXT","w");
	phis2=0.;
	xnt=96.6;
	w=2.;
	phasedeg=0.;
	sigrin=.0001;
	siggl=0.;
	srn=0.;
	ra=21000.;
	whic=-1.;
	ts=.01;
	tf=10.;
	phis1=w*w*xnt*xnt/tf;
	qperfect=0;
	vc=9000.;
	xnp=3.;
	xnclim=9999999.;
	apn=4;
	tau=.5;
	hedeg=0.;
	vm=3000.;
	phase=phasedeg/57.3;
	order=5;
	x=w*ts;
	tgo=tf;
	t=0.;
	s=0.;
	y=0.;
	yd=-xnt/w-vm*hedeg/57.3;
	ytdd=xnt*sin(w*t);
	ytddd=xnt*w*cos(w*t);
	xnc=0.;
	xnl=0.;
	h=.001;
	hp=.001;
	ts2=ts*ts;
	ts3=ts2*ts;
	ts4=ts3*ts;
	ts5=ts4*ts;
	ts6=ts5*ts;
	ts7=ts6*ts;
	wh=whic;
	if(qperfect==1){
		yh=y;
		ydh=yd;
		ytddh=ytdd;
		ytdddh=ytddd;
		wh=w;
	}
	else{
		yh=0.;
		ydh=0.;
		ytddh=0.;
		ytdddh=0.;
	}
	for (i=1; i<=order; i=i+1){
		for (j=1; j<=order; j=j+1){
			phi[i][j]=0.;
			p[i][j]=0.;
			q[i][j]=0.;
			idn[i][j]=0.;
		}
	}
 	idn[1][1]=1.;
	idn[2][2]=1.;
	idn[3][3]=1.;
	idn[4][4]=1.;
	idn[5][5]=1.;
	rtm=vc*tf;
	signoise=sigrin;
	ynoise=signoise*rtm;
	p[1][1]=ynoise*ynoise;
	p[2][2]=(vm*20./57.3)*(vm*20./57.3);
	p[3][3]=xnt*xnt;
	p[4][4]=(w*xnt)*(w*xnt);
	p[5][5]=w*w;
	for (i=1; i<=order; i=i+1){
		hmat[1][i]=0.;
		ht[i][1]=0.;
	}
	hmat[1][1]=1.;
	ht[1][1]=1.;
L10:
	if(t>(tf-.0001))goto L999;
 	yold=y;
	ydold=yd;
	xnlold=xnl;
	step=1;
	goto L200;
L66:
	step=2;
 	y=y+h*yd;
 	yd=yd+h*ydd;
	xnl=xnl+h*xnld;
	t=t+h;
	goto L200;
L55:
 	y=.5*(yold+y+h*yd);
 	yd=.5*(ydold+yd+h*ydd);
	xnl=.5*(xnlold+xnl+h*xnld);
	s=s+h;
	if(s<(ts-.00001))goto L10;
	s=0.;
	ytdd=xnt*sin(w*t);
	ytddd=xnt*w*cos(w*t);
	phi[1][1]=1.;
	phi[1][2]=ts;
	phi[2][2]=1.;
	phi[2][3]=ts;
	phi[3][3]=1.;
	phi[3][4]=ts;
	phi[4][3]=-wh*wh*ts;
	phi[4][4]=1.;
	phi[4][5]=-2.*wh*ytddh*ts;
	phi[5][5]=1.;
	q[3][3]=phis1*ts*ts*ts/3.;
	q[3][4]=phis1*ts*ts/2.;
	q[4][3]=q[3][4];
	q[4][4]=4.*wh*wh*ytddh*ytddh*phis2*ts*ts*ts/3.+phis1*ts;
	q[4][5]=-wh*ytddh*ts*ts*phis2;
	q[5][4]=q[4][5];
	q[5][5]=phis2*ts;
	for (i=1; i<=order; i=i+1){
		for (j=1; j<=order; j=j+1){
			phit[j][i]=phi[i][j];
		}
	}
	for (i=1; i<=order; i=i+1){
		for (j=1; j<=order; j=j+1){
			phip[i][j]=0.;
			for (k=1; k<=order; k=k+1){
				phip[i][j]=phip[i][j]+phi[i][k]*p[k][j];
			}
		}
	}
	for (i=1; i<=order; i=i+1){
		for (j=1; j<=order; j=j+1){
			phipphit[i][j]=0.;
			for (k=1; k<=order; k=k+1){
				phipphit[i][j]=phipphit[i][j]+phip[i][k]*phit[k][j];
			}
		}
	}
	for (i=1; i<=order; i=i+1){
		for (j=1; j<=order; j=j+1){
			m[i][j]=phipphit[i][j]+q[i][j];
		}
	}
	for (i=1; i<=1; i=i+1){
		for (j=1; j<=order; j=j+1){
			hm[i][j]=0.;
			for (k=1; k<=order; k=k+1){
				hm[i][j]=hm[i][j]+hmat[i][k]*m[k][j];
			}
		}
	}
 	for (i=1; i<=1; i=i+1){
		for (j=1; j<=1; j=j+1){
			hmht[i][j]=0.;
			for (k=1; k<=order; k=k+1){
				hmht[i][j]=hmht[i][j]+hm[i][k]*ht[k][j];
			}
		}
	}
 	rtm=vc*tgo;
	signoise=sigrin;
	ynoise=signoise*rtm;
	rmat[1][1]=ynoise*ynoise;
	hmhtr[1][1]=hmht[1][1]+rmat[1][1];
	hmhtrinv[1][1]=1./hmhtr[1][1];
	for (i=1; i<=order; i=i+1){
		for (j=1; j<=1; j=j+1){
			mht[i][j]=0.;
			for (k=1; k<=order; k=k+1){
				mht[i][j]=mht[i][j]+m[i][k]*ht[k][j];
			}
		}
	}
	for (i=1; i<=order; i=i+1){
		gain[i][1]=mht[i][1]*hmhtrinv[1][1];
	}
	for (i=1; i<=order; i=i+1){
		for (j=1; j<=order; j=j+1){
			kh[i][j]=0.;
			for (k=1; k<=1; k=k+1){
				kh[i][j]=kh[i][j]+gain[i][k]*hmat[k][j];
			}
		}
	}
	for (i=1; i<=order; i=i+1){
		for (j=1; j<=order; j=j+1){
			ikh[i][j]=idn[i][j]-kh[i][j];
		}
	}
 	for (i=1; i<=order; i=i+1){
		for (j=1; j<=order; j=j+1){
			p[i][j]=0.;
			for (k=1; k<=order; k=k+1){
				p[i][j]=p[i][j]+ikh[i][k]*m[k][j];
			}
		}
	}
 	rtm=vc*tgo;
	xlam=y/rtm;
	gauss(signoise,&xnoise);
	xlams=xlam+xnoise;
	project(ts,yh,ydh,ytddh,ytdddh,&yb,&ydb,&ytddb,&ytdddb,hp,xnl,wh);
	res=rtm*xlams-yb;
	yh=yb+gain[1][1]*res;
	ydh=ydb+gain[2][1]*res;
	ytddh=ytddb+gain[3][1]*res;
	ytdddh=ytdddb+gain[4][1]*res;
	wh=wh+gain[5][1]*res;
	if(apn==0){
		xnc=xnp*(yh+ydh*tgo)/(tgo*tgo);
	}
	else if(apn==1){
		xnc=xnp*(yh+ydh*tgo+.5*ytddh*tgo*tgo)/(tgo*tgo);
	}
	else if(apn==2){
		xs=tgo/tau;
		top=6.*xs*xs*(exp(-xs)-1.+xs);
		bot1=2*xs*xs*xs+3.+6.*xs-6.*xs*xs;
		bot2=-12.*xs*exp(-xs)-3.*exp(-2.*xs);
		xnpp=top/(.0001+bot1+bot2);
		c1=xnpp/(tgo*tgo);
		c2=xnpp/tgo;
		c3=.5*xnpp;
		c4=-xnpp*(exp(-xs)+xs-1.)/(xs*xs);
		xnc=c1*yh+c2*ydh+c3*ytddh+c4*xnl;
	}
	else if(apn==3){
		xp=wh*tgo;
		xnc=xnp*(yh+ydh*tgo)/(tgo*tgo)+xnp*ytddh*(1.-cos(xp))/(xp*xp)+
     			xnp*ytdddh*(xp-sin(xp))/(xp*xp*wh);
     	}
	else{
		xs=tgo/tau;
		top=6.*xs*xs*(exp(-xs)-1.+xs);
		bot1=2*xs*xs*xs+3.+6.*xs-6.*xs*xs;
		bot2=-12.*xs*exp(-xs)-3.*exp(-2.*xs);
		xnpp=top/(.0001+bot1+bot2);
		c1=xnpp/(tgo*tgo);
		c2=xnpp/tgo;
		c3=xnpp*(1.-cos(wh*tgo))/(wh*wh*tgo*tgo);
		c4=-xnpp*(exp(-xs)+xs-1.)/(xs*xs);
		c5=xnpp*(wh*tgo-sin(wh*tgo))/(wh*wh*wh*tgo*tgo);
		xnc=c1*yh+c2*ydh+c3*ytddh+c4*xnl+c5*ytdddh;
	}
	if(xnc>xnclim)xnc=xnclim;
	if(xnc<-xnclim)xnc=-xnclim;
	errytdd=ytdd-ytddh;
	errytddg=errytdd/32.2;
	errytddd=ytddd-ytdddh;
	errytdddg=errytddd/32.2;
	sp55=sqrt(p[5][5]);
	sp55p=-sp55;
	errw=w-wh;
	sp44=sqrt(p[4][4]);
	sp44p=-sp44;
	sp33=sqrt(p[3][3]);
	sp33p=-sp33;
	sp33g=sp33/32.2;
	sp33pg=sp33p/32.2;
	sp44g=sp44/32.2;
	sp44pg=sp44p/32.2;
	ytddg=ytdd/32.2;
	ytddhg=ytddh/32.2;
	ytdddg=ytddd/32.2;
	ytdddhg=ytdddh/32.2;
	xncg=xnc/32.2;
	printf("%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n"
		,t,ytddg,ytddhg,ytdddg,ytdddhg,w,wh);
	fprintf(fptr1,"%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n"
		,t,ytddg,ytddhg,ytdddg,ytdddhg,w,wh);
	fprintf(fptr2,"%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n"
		,t,errytdd,sp33,-sp33,errytddd,sp44,-sp44,errw,sp55,-sp55);
     	goto L10;
L200:
	ytdd=xnt*sin(w*t);
 	tgo=tf-t+.00001;
	xnld=(xnc-xnl)/tau;
	ydd=ytdd-xnl;
	if (step<2)
		goto L66;
	else
		goto L55;
L999:
 	fclose (fptr1);
 	fclose (fptr2);
 	printf("%6.3f\n",y);
 	return 0;
}
	
	void gauss(double signoise,double *xnoise)
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
	*xnoise=1.414*temp*signoise;
}
	
	
void project(double ts,double yph,double ydph,double ytddph,double ytdddph,
	double *yb,double *ydb,double *ytddb,double *ytdddb,double hp,double xnl,
	double wph)
{
	double t,y,yd,ydd,ytddd,w,h,ytdddd,ytdd;
	t=0.;
	y=yph;
	yd=ydph;
	ytdd=ytddph;
	ytddd=ytdddph;
	w=wph;
	h=hp;
L10:
	if(t>(ts-.0001))goto L999;
	ytdddd=-w*w*ytdd;
	ytddd=ytddd+h*ytdddd;
	ytdd=ytdd+h*ytddd;
	ydd=ytdd-xnl;
	yd=yd+h*ydd;
	y=y+h*yd;
	t=t+h;
	goto L10;
L999:
 	*yb=y;
 	*ydb=yd;
	*ytddb=ytdd;
	*ytdddb=ytddd;
}
	
