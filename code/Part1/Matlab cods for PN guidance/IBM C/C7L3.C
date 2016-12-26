#include <stdio.h>
#include <stdlib.h>
#include <math.h>
main()
{
	int step;
	float xnt,xnp,tf,ts,beta,signoise,vc,t,s,tp;
	float x1,x2,x3,x5,y1old,y2old,y3old,y4old,y5old,h;
	float x1old,x2old,x3old,x5old;
	float gfilter,hfilter,x1d,x2d,x3d,x5d,tgo;
	float temp1,temp2,y1new,y2new,y3new,y4new,y5new,xmnoise,xmnt;
	FILE *fptr;
	fptr=fopen("DATFIL.TXT","w");
	xnt=96.6;
	xnp=3.;
	tf=10.;
	ts=.1;
	beta=.8;
	signoise=.001;
	vc=4000.;
	t=0.;
	s=0.;
	tp=t+.00001;
	x1=0;
	x2=0;
	x3=1;
	x5=0.;
	y1old=0.;
	y2old=0.;
	y3old=0.;
	y4old=0.;
	y5old=0.;
	h=.01;
	gfilter=1.-beta*beta;
	hfilter=(1.-beta)*(1.-beta);
L10:	if(tp>(tf-.00001))goto L999;
	s=s+h;
	x1old=x1;
	x2old=x2;
	x3old=x3;
	x5old=x5;
	step=1;
	goto L200;
L66:	step=2;
	x1=x1+h*x1d;
	x2=x2+h*x2d;
	x3=x3+h*x3d;
	x5=x5+h*x5d;
	tp=tp+h;
	goto L200;
L55:
	x1=(x1old+x1)/2+.5*h*x1d;
	x2=(x2old+x2)/2+.5*h*x2d;
	x3=(x3old+x3)/2+.5*h*x3d;
	x5=(x5old+x5)/2+.5*h*x5d;
	if(s<(ts-.0001))goto L10;
	s=0.;
	temp1=(x5-y1old)*xnp*vc;
	temp2=hfilter*(y2old+temp1)/ts+gfilter*y3old;
	y1new=x5;
	y2new=temp1+y2old+ts*(y3old-temp2);
	y3new=y3old-temp2;
	y4new=y4old+temp2;
	y5new=y5old+temp2*temp2;
	y1old=y1new;
	y2old=y2new;
	y3old=y3new;
	y4old=y4new;
	y5old=y5new;
	xmnoise=signoise*sqrt(y5new);
	xmnt=xnt*x1;
	printf("%6.3f %6.3f %6.3f\n",tp,xmnt,xmnoise);
	fprintf(fptr,"%6.3f %6.3f %6.3f\n",tp,xmnt,xmnoise);
	goto L10;
L200:
 	tgo=tp+.00001;
	x1d=x2;
	x2d=x3+y4old/(vc*tgo);
	x3d=(y4old)/(vc*tgo*tgo);
	x5d=-x2;
	if (step<2)
		goto L66;
	else
		goto L55;
L999:
	fclose (fptr);
	return 0;
}
