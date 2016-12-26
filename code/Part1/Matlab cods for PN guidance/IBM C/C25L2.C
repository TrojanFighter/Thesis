#include <stdio.h>
#include <stdlib.h>
#include <math.h>
main()
{
	int step,apn;
	float xntg,hedeg,xnp,rm1ic,rm2ic,rt1ic,rt2ic,vm,vt,xnclimg;
	float xlamfdeg,xnclim,xlamf,xnt,rm1,rm2,rt1,rt2,beta;
	float vt1,vt2,he,t,s,rtm1,rtm2,rtm,xlam,xlead,thet,vm1,vm2;
	float vtm1,vtm2,vc,h,betaold,rt1old,rt2old,rm1old,rm2old;
	float vm1old,vm2old,rt1k,rt2k,rm1k,rm2k,appang,xncg,xlamd,tgo,xnc,xnt1,xnt2;
	float xntplos,am1,am2,betad;
	FILE *fptr;
	fptr=fopen("DATFIL.TXT","w");
	xntg=0.;
	hedeg=0.;
	xnp=3.;
	rm1ic=0.;
	rm2ic=10000.;
	rt1ic=30000.;
	rt2ic=0.;
	vm=3000.;
	vt=0.;
	xnclimg=9999999.;
	apn=1;
	xlamfdeg=-90.;
	h=.0001;
	xnclim=32.2*xnclimg;
	xlamf=xlamfdeg/57.3;
	xnt=32.2*xntg;
	rm1=rm1ic;
	rm2=rm2ic;
	rt1=rt1ic;
	rt2=rt2ic;
	beta=0.;
	vt1=-vt*cos(beta);
	vt2=vt*sin(beta);
	he=hedeg/57.3;
	t=0.;
	s=0.;
	rtm1=rt1-rm1;
	rtm2=rt2-rm2;
	rtm=sqrt(rtm1*rtm1+rtm2*rtm2);
	xlam=atan2(rtm2,rtm1);
	xlead=asin(vt*sin(beta+xlam)/vm);
	thet=xlam+xlead;
	vm1=vm*cos(thet+he);
	vm2=vm*sin(thet+he);
	vtm1=vt1-vm1;
	vtm2=vt2-vm2;
	vc=-(rtm1*vtm1+rtm2*vtm2)/rtm;
L10:
	if(vc<0.)goto L999;
 	if(rtm<1000.)
		h=.00001;
	else
		h=.0001;
 	betaold=beta;
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
 	beta=beta+h*betad;
 	rt1=rt1+h*vt1;
	rt2=rt2+h*vt2;
	rm1=rm1+h*vm1;
	rm2=rm2+h*vm2;
	vm1=vm1+h*am1;
	vm2=vm2+h*am2;
	t=t+h;
	goto L200;
L55:
 	beta=.5*(betaold+beta+h*betad);
 	rt1=.5*(rt1old+rt1+h*vt1);
	rt2=.5*(rt2old+rt2+h*vt2);
	rm1=.5*(rm1old+rm1+h*vm1);
	rm2=.5*(rm2old+rm2+h*vm2);
	vm1=.5*(vm1old+vm1+h*am1);
	vm2=.5*(vm2old+vm2+h*am2);
	s=s+h;
	if(s<.09999)goto L10;
	s=0.;
	rt1k=rt1/1000.;
	rt2k=rt2/1000.;
	rm1k=rm1/1000.;
	rm2k=rm2/1000.;
	appang=180.+xlam*57.3;
	xncg=xnc/32.2;
	printf("%10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n",t,rm1k,rm2k,rt1k,rt2k,xnc/32.2,xlam*57.3);
	fprintf(fptr,"%10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n",t,rm1k,rm2k,rt1k,rt2k,xnc/32.2,xlam*57.3);
	goto L10;
L200:
 	rtm1=rt1-rm1;
	rtm2=rt2-rm2;
	rtm=sqrt(rtm1*rtm1+rtm2*rtm2);
	vtm1=vt1-vm1;
	vtm2=vt2-vm2;
	vc=-(rtm1*vtm1+rtm2*vtm2)/rtm;
	xlam=atan2(rtm2,rtm1);
	xlamd=(rtm1*vtm2-rtm2*vtm1)/(rtm*rtm);
	tgo=rtm/vc;
	if(apn==0)
		xnc=xnp*vc*xlamd;
	else
		xnt1=xnt*sin(beta);
		xnt2=xnt*cos(beta);
		xntplos=-xnt1*sin(xlam)+xnt2*cos(xlam);
		xnc=4.*vc*xlamd+xntplos+2.*vc*(xlam-xlamf)/tgo;
	if(xnc>xnclim)xnc=xnclim;
	if(xnc<-xnclim)xnc=-xnclim;
	am1=-xnc*sin(xlam);
	am2=xnc*cos(xlam);
	vt1=-vt*cos(beta);
	vt2=vt*sin(beta);
	if(vt==0.)
		betad=0.;
	else
		betad=xnt/vt;
	if(step<=1)
		goto L66;
	else
		goto L55;
L999:
 	printf("%10.4f\n",rtm);
	fclose (fptr);
	return 0;
}

	
