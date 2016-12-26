#include <stdio.h>
#include <stdlib.h>
#include <math.h>
main()
{	
	int step;
	double xt,yt,tf,h,ts,xnc,xm,ym,xmd,ymd,s,t,xmold,ymold,xmdold,ymdold;
	double tgo,rtm1,rtm2,rtm,vtm1,vtm2,phi,xmdd,ymdd;
	FILE *fptr1;
	fptr1=fopen("DATFIL.TXT","w");
	xt=-4000.;
	yt=5000.;
	tf=100.;
	h=.01;
	ts=.1;
	xnc=12.;
	xm=0.;
	ym=0.;
	xmd=0.;
	ymd=0.;
	s=0.;
	t=0.;
	phi=0.;
L10:	
	if(t>tf)goto L999;
	xmold=xm;
	ymold=ym;
	xmdold=xmd;
	ymdold=ymd;
	step=1;
	goto L200;
L66:
	step=2;
	xm=xm+h*xmd;
	ym=ym+h*ymd;
	xmd=xmd+h*xmdd;
	ymd=ymd+h*ymdd;
	t=t+h;
	goto L200;
L55:
	xm=(xmold+xm)/2+.5*h*xmd;
	ym=(ymold+ym)/2+.5*h*ymd;
	xmd=(xmdold+xmd)/2+.5*h*xmdd;
	ymd=(ymdold+ymd)/2+.5*h*ymdd;
	s=s+h;
	if(s<(ts-.0001))goto L10;
 	s=0.;
	tgo=tf-t+.0001;
	rtm1=xt-xm;
	rtm2=yt-ym;
	rtm=sqrt(rtm1*rtm1+rtm2*rtm2);
	vtm1=-xmd;
	vtm2=-ymd;
	phi=atan2(rtm2+vtm2*tgo,rtm1+vtm1*tgo);
	printf("%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n"
		,t,xm,ym,xt,yt,phi*57.3);
	fprintf(fptr1,"%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n"
		,t,xm,ym,xt,yt,phi*57.3);
 	goto L10;
L200:
 	xmdd=xnc*cos(phi);
	ymdd=xnc*sin(phi);
	if (step<2)
		goto L66;
	else
		goto L55;
L999:
 	rtm1=xt-xm;
	rtm2=yt-ym;
	rtm=sqrt(rtm1*rtm1+rtm2*rtm2);
 	printf("%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n"
		,t,xm,ym,xt,yt,phi*57.3);
	fprintf(fptr1,"%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n"
		,t,xm,ym,xt,yt,phi*57.3);
	printf("%6.3f\n",rtm);
	fclose (fptr1);
	return 0;
}
	
