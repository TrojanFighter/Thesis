#include <stdio.h>
#include <stdlib.h>
#include <math.h>
void gauss(float signoise,float *xlamnoise);
main()
{
	float vc,xnt,yic,vm,hedeg,beta,xnp,signoise,tf,ts;
	float y,yd,ydic,t,h,s,gfilter,hfilter,xlamh,xlamdh,xnc;
	float yold,ydold,ydd,res,xlamnoise;
	float tgo,rtm,xlam,xlamd;
	int noise,step;
	FILE *fptr;
	fptr=fopen("DATFIL.TXT","w");
	vc=4000.;
	xnt=96.6;
	yic=0.;
	vm=3000.;
	hedeg=0.;
	beta=.3;
	xnp=3.;
	signoise=.001;
	tf=10.;
	ts=.1;
	noise=1;
	y=yic;
	yd=-vm*hedeg/57.3;
	ydic=yd;
	t=0.;
	h=.01;
	s=0.;
	gfilter=1.-beta*beta;
	hfilter=(1.-beta)*(1.-beta);
	xlamh=0.;
	xlamdh=0.;
	xnc=0.;
L10:
	if(t>(tf-.0001))goto L999;
 	yold=y;
	ydold=yd;
	step=1;
	goto L200;
L66:
	step=2;
 	y=y+h*yd;
 	yd=yd+h*ydd;
	t=t+h;
	goto L200;
L55:
 	y=.5*(yold+y+h*yd);
 	yd=.5*(ydold+yd+h*ydd);
	s=s+h;
	if(s<(ts-.0001))goto L10;
	s=0.;
	if(noise==1)
		gauss(signoise,&xlamnoise);
	else
		xlamnoise=0.;
	res=xlam-(xlamh+ts*xlamdh)+xlamnoise;
	xlamh=gfilter*res+xlamh+ts*xlamdh;
	xlamdh=hfilter*res/ts+xlamdh;
	xnc=xnp*vc*xlamdh;
	y=fabs(y);
	printf("%6.3f %6.3f %6.3f %6.3f %6.3f\n", t,y,xnc,xlamd,xlamdh);
	fprintf(fptr,"%6.3f %6.3f %6.3f %6.3f %6.3f\n", t,y,xnc,xlamd,xlamdh);
	goto L10;
L200:
 	tgo=tf-t+.00001;
	rtm=vc*tgo;
	xlam=y/(vc*tgo);
	xlamd=(rtm*yd+y*vc)/(rtm*rtm);
	ydd=xnt-xnc;
	if (step<2)
		goto L66;
	else
		goto L55;
L999:
 	y=fabs(y);
 	printf("%6.3f %6.3f\n", t,rtm);
	fclose (fptr);
	return 0;
}

	void gauss(float signoise,float *xlamnoise)
{
	float sum,x,y,temp;
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
