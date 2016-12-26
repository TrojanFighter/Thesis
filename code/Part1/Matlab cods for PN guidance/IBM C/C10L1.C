#include <stdio.h>
#include <stdlib.h>
#include <math.h>
void gauss(float signoise,float *thettnoise);
main()
{
	float vm,vt,xnt,rm1ic,rm2ic,rt1ic,rt2ic,hedeg,xnp,beta,ts,signoise;
	int noise,step;
	float rt1,rt2,rm1,rm2,betat,vt1,vt2,he,gfilter,hfilter,h,thetm;
	float xlamh,xlamdh,xnc,t,s,rtm1,rtm2,rtm,xlam,xlead,xlamd;
	float thet,vm1,vm2,vtm1,vtm2,vc,betatold,rt1old,rt2old,rm1old,rm2old;
	float vm1old,vm2old,am1,am2,betatd,thettm,thett,thetmm,thettnoise;
	float rt1m,rt2m,rm1m,rm2m,xlamm,res,rt1km,rt2km,rm1km,rm2km,rt,rm;
	FILE *fptr;
	fptr=fopen("DATFIL.TXT","w");
	vm=3000.;
	vt=1000.;
	xnt=0.;
	rm1ic=0.;
	rm2ic=10000.;
	rt1ic=40000.;
	rt2ic=10000.;
	hedeg=0.;
	xnp=3.;
	beta=.8;
	ts=.1;
	signoise=.001;
	noise=1;
	rt1=rt1ic;
	rt2=rt2ic;
	rm1=rm1ic;
	rm2=rm2ic;
	betat=0.;
	vt1=-vt*cos(betat);
	vt2=vt*sin(betat);
	he=hedeg/57.3;
	gfilter=1.-beta*beta;
	hfilter=(1.-beta)*(1.-beta);
	xlamh=0.;
	xlamdh=0.;
	xnc=0.;
	t=0.;
	s=0.;
	rtm1=rt1-rm1;
	rtm2=rt2-rm2;
	rtm=sqrt(rtm1*rtm1+rtm2*rtm2);
	xlam=atan2(rtm2,rtm1);
	xlead=asin(vt*sin(betat+xlam)/vm);
	thet=xlam+xlead;
	vm1=vm*cos(thet+he);
	vm2=vm*sin(thet+he);
	vtm1=vt1-vm1;
	vtm2=vt2-vm2;
	vc=-(rtm1*vtm1+rtm2*vtm2)/rtm;
L10:
	if(vc<0.)goto L999;
 	if(rtm<1000.)
		h=.0002;
	else
		h=.01;
 	betatold=betat;
	rt1old=rt1;
	rt2old=rt2;
	rm1old=rm1;
	rm2old=rm2;
	vm1old=vm1;
	vm2old=vm2;
	step=1;
	goto L200;
L66:
	step=2;
 	betat=betat+h*betatd;
 	rt1=rt1+h*vt1;
	rt2=rt2+h*vt2;
	rm1=rm1+h*vm1;
	rm2=rm2+h*vm2;
	vm1=vm1+h*am1;
	vm2=vm2+h*am2;
	t=t+h;
	goto L200;
L55:
 	betat=.5*(betatold+betat+h*betatd);
 	rt1=.5*(rt1old+rt1+h*vt1);
	rt2=.5*(rt2old+rt2+h*vt2);
	rm1=.5*(rm1old+rm1+h*vm1);
	rm2=.5*(rm2old+rm2+h*vm2);
	vm1=.5*(vm1old+vm1+h*am1);
	vm2=.5*(vm2old+vm2+h*am2);
	s=s+h;
	if(s<(ts-.0001))goto L10;
	s=0.;
	if(noise==1)
		gauss(signoise,&thettnoise);
	else
		thettnoise=0.;
	thettm=thett+thettnoise;
	thetmm=thetm;
	rt1m=rt*cos(thettm);
	rt2m=rt*sin(thettm);
	rm1m=rm*cos(thetmm);
	rm2m=rm*sin(thetmm);
	xlamm=atan2(rt2m-rm2m,rt1m-rm1m);
	res=xlamm-(xlamh+ts*xlamdh);
	xlamh=gfilter*res+xlamh+ts*xlamdh;
	xlamdh=hfilter*res/ts+xlamdh;
	xnc=xnp*vc*xlamdh;
	rt1km=rt1/3280.;
	rt2km=rt2/3280.;
	rm1km=rm1/3280.;
	rm2km=rm2/3280.;
	printf("%6.3f %8.6f %8.6f %6.3f %6.3f\n", t,xlamd,xlamdh,xnc,rtm);
	fprintf(fptr,"%6.3f %8.6f %8.6f %6.3f %6.3f\n", t,xlamd,xlamdh,xnc,rtm);
	goto L10;
L200:
 	thett=atan2(rt2,rt1);
	thetm=atan2(rm2,rm1);
	rt=sqrt(rt1*rt1+rt2*rt2);
	rm=sqrt(rm1*rm1+rm2*rm2);
 	rtm1=rt1-rm1;
	rtm2=rt2-rm2;
	rtm=sqrt(rtm1*rtm1+rtm2*rtm2);
	vtm1=vt1-vm1;
	vtm2=vt2-vm2;
	vc=-(rtm1*vtm1+rtm2*vtm2)/rtm;
	xlam=atan2(rtm2,rtm1);
	xlamd=(rtm1*vtm2-rtm2*vtm1)/(rtm*rtm);
	xnc=xnp*vc*xlamdh;
	am1=-xnc*sin(xlam);
	am2=xnc*cos(xlam);
	vt1=-vt*cos(betat);
	vt2=vt*sin(betat);
	betatd=xnt/vt;
	if (step<2)
		goto L66;
	else
		goto L55;
L999:
 	printf("%6.3f %6.3f %6.3f %6.3f %6.3f\n", t,xlamd,xlamdh,xnc,rtm);
	fclose (fptr);
	return 0;
}
	
	void gauss(float signoise,float *thettnoise)
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
	*thettnoise=1.414*temp*signoise;
}
