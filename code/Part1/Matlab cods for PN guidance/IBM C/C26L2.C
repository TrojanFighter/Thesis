#include <stdio.h>
#include <stdlib.h>
#include <math.h>
void gauss(double signoise,double *xlamnoise);
main()
{	
	double phi[5][5],p[5][5],m[5][5],phip[5][5],phipphit[5][5],gain[5][2];
	double q[5][5],hmat[2][5],hm[2][5],mht[5][2];
	double phit[5][5];
	double hmht[2][2],ht[5][2],kh[5][5],idn[5][5],ikh[5][5];
	int apn,order,mvr,qperfect,i,j,k,step;
	double tau,vc,w,wreal,wh,xnt,xntreal,ts,yic,vm,hedeg,hedegfil,xnp;
	double sigrin,tf,phase,x,y,yd,phis,rtm,signoise,sigpos;
	double sign2,t,h,s,xnc,xnl,xlam,yold,ydold,xnlold,tgo,xlamnoise;
	double ystar,res,ytddhnew,ytddg,ytddhg,erry,sp11,sp11p,erryd,sp22,sp22p,errytddg,sp33g;
	double errytdddg,sp44g,xnld,ydd,ytdd,ytddd,yh,ydh,ytddh,ytdddh,hmhtr,hmhtrinv;
	double xs,top,bot1,bot2,xnpp,c1,c2,c3,c4,xp,c5;
	FILE *fptr1;
	FILE *fptr2;
	fptr1=fopen("DATFIL.TXT","w");
	fptr2=fopen("COVFIL.TXT","w");
	tau=.5;
	apn=0;
	order=4;
	mvr=1;
	vc=9000.;
	w=2.;
	wreal=2.;
	wh=w;
	xnt=96.6;
	xntreal=96.6;
	ts=.01;
	yic=0.;
	vm=3000.;
	hedeg=0.;
	hedegfil=20.;
	xnp=3.;
	sigrin=.001;
	tf=10.;
	qperfect=0;
	phase=0./57.3;
	x=wh*ts;
	y=yic;
	yd=-vm*hedeg/57.3;
	phis=wh*wh*xnt*xnt/tf;
	rtm=vc*tf;
	signoise=sigrin;
	sigpos=rtm*signoise;
	sign2=sigpos*sigpos;
	for (i=1; i<=order; i=i+1){
		for (j=1; j<=order; j=j+1){
			phi[i][j]=0.;
			p[i][j]=0.;
			q[i][j]=0.;
			idn[i][j]=0.;
		}
	}
	phi[1][1]=1;
	phi[1][2]=ts;
	phi[1][3]=(1-cos(x))/(wh*wh);
	phi[1][4]=(x-sin(x))/(wh*wh*wh);
	phi[2][2]=1;
	phi[2][3]=sin(x)/wh;
	phi[2][4]=(1-cos(x))/(wh*wh);
	phi[3][3]=cos(x);
	phi[3][4]=sin(x)/wh;
	phi[4][3]=-wh*sin(x);
	phi[4][4]=cos(x);
	q[1][1]=phis*(.333*x*x*x-2*sin(x)+2*x*cos(x)+.5*x-.25*sin(2*x))/(wh*wh*wh*wh*wh);
	q[1][2]=phis*(.5*x*x-x*sin(x)+.5*sin(x)*sin(x))/(wh*wh*wh*wh);
	q[2][1]=q[1][2];
	q[1][3]=phis*(sin(x)-x*cos(x)-.5*x+.25*sin(2*x))/(wh*wh*wh);
	q[3][1]=q[1][3];
	q[1][4]=phis*(cos(x)+x*sin(x)-.5*sin(x)*sin(x)-1)/(wh*wh);
	q[4][1]=q[1][4];
	q[2][2]=phis*(1.5*x-2*sin(x)+.25*sin(2*x))/(wh*wh*wh);
	q[2][3]=phis*(1-cos(x)-.5*sin(x)*sin(x))/(wh*wh);
	q[3][2]=q[2][3];
	q[2][4]=phis*(sin(x)-.5*x-.25*sin(2*x))/wh;
	q[4][2]=q[2][4];
	q[3][3]=phis*(.5*x-.25*sin(2*x))/wh;
	q[3][4]=.5*phis*sin(x)*sin(x);
	q[4][3]=q[3][4];
	q[4][4]=wh*phis*(.5*x+.25*sin(2*x));
 	idn[1][1]=1.;
	idn[2][2]=1.;
	idn[3][3]=1.;
	idn[4][4]=1.;
	p[1][1]=sign2;
	p[2][2]=(vm*hedegfil/57.3)*(vm*hedegfil/57.3);
	p[3][3]=xnt*xnt;
	p[4][4]=wh*wh*xnt*xnt;
	for (i=1; i<=order; i=i+1){
		hmat[1][i]=0.;
		ht[i][1]=0.;
	}
	hmat[1][1]=1.;
	ht[1][1]=1.;
	for (i=1; i<=order; i=i+1){
		for (j=1; j<=order; j=j+1){
			phit[j][i]=phi[i][j];
		}
	}
	t=0.;
	h=.001;
	s=0.;
	xnc=0.;
	xnl=0.;
	xlam=y/rtm;
	if(mvr==0){
		ytdd=xntreal;
		ytddd=0.;
	}
	else{
		ytdd=xntreal*sin(wreal*t);
		ytddd=xntreal*wreal*cos(wreal*t);
	}
	if(qperfect==1){
		yh=y;
		ydh=yd;
		ytddh=ytdd;
		ytdddh=ytddd;
	}
	else{
		yh=0.;
		ydh=0.;
		ytddh=0.;
		ytdddh=0.;
	}
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
	tgo=tf-t+.000001;
	rtm=vc*tgo;
	signoise=sigrin;
	sigpos=rtm*signoise;
	sign2=sigpos*sigpos;
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
	hmhtr=hmht[1][1]+sign2;
	hmhtrinv=1./hmhtr;
	for (i=1; i<=order; i=i+1){
		for (j=1; j<=1; j=j+1){
			mht[i][j]=0.;
			for (k=1; k<=order; k=k+1){
				mht[i][j]=mht[i][j]+m[i][k]*ht[k][j];
			}
		}
	}
	for (i=1; i<=order; i=i+1){
		gain[i][1]=mht[i][1]*hmhtrinv;
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
 	if(mvr==0){
		ytdd=xntreal;
		ytddd=0.;
	}
	else{
		ytdd=xntreal*sin(wreal*t);
		ytddd=xntreal*wreal*cos(wreal*t);
	}
	gauss(signoise,&xlamnoise);
	ystar=rtm*(xlam+xlamnoise);
	res=ystar-yh-ts*ydh-(1-cos(x))*ytddh/(wh*wh)-(x-sin(x))*ytdddh/(wh*wh*wh)+.5*ts*ts*xnl;
     	yh=yh+ts*ydh+(1-cos(x))*ytddh/(wh*wh)+(x-sin(x))*ytdddh/(wh*wh*wh)+gain[1][1]*res-.5*ts*ts*xnl;
     	ydh=ydh+sin(x)*ytddh/wh+(1-cos(x))*ytdddh/(wh*wh)+gain[2][1]*res-ts*xnl;
     	ytddhnew=cos(x)*ytddh+sin(x)*ytdddh/wh+gain[3][1]*res;
	ytdddh=-wh*sin(x)*ytddh+cos(x)*ytdddh+gain[4][1]*res;
	ytddh=ytddhnew;
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
		xnc=xnp*(yh+ydh*tgo)/(tgo*tgo)+xnp*ytddh*(1.-cos(xp))/xp*xp+xnp*ytdddh*(xp-sin(xp))/(xp*xp*wh);
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
	ytddg=ytdd/32.2;
	ytddhg=ytddh/32.2;
	erry=y-yh;
	sp11=sqrt(p[1][1]);
	sp11p=-sp11;
	erryd=yd-ydh;
	sp22=sqrt(p[2][2]);
	sp22p=-sp22;
	errytddg=(ytdd-ytddh)/32.2;
	sp33g=sqrt(p[3][3])/32.2;
	errytdddg=(ytddd-ytdddh)/32.2;
	sp44g=sqrt(p[4][4])/32.2;
	printf("%6.3f %6.3f %6.3f %6.3f %6.3f\n"
			,t,ytdd/32.2,ytddh/32.2,ytddd/32.2,ytdddh/32.2);
	fprintf(fptr1,"%6.3f %6.3f %6.3f %6.3f %6.3f\n"
			,t,ytdd/32.2,ytddh/32.2,ytddd/32.2,ytdddh/32.2);
	fprintf(fptr2,"%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n"
			,t,errytddg,sp33g,-sp33g,errytdddg,sp44g,-sp44g);
	goto L10;
L200:
 	tgo=tf-t+.000001;
	rtm=vc*tgo;
	xlam=y/(vc*tgo);
	if(mvr==0)
		ytdd=xntreal;
	else
		ytdd=xntreal*sin(wreal*t);
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
