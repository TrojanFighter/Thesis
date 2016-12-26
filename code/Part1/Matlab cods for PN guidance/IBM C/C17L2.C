#include <stdio.h>
#include <stdlib.h>
#include <math.h>
void initial(float rt2des,float rt1,float rt2,float vt1,float vt2,float betest,
float *rt1f,float *rt2f,float *tfdes);
main()
{
	int step,apn;
	float xnp,rt1,rt2,rm1,rm2,vt,rt2des,gamtdeg,beta,betest;
	float xnclimg,xnclim,vt1,vt2,rt1f,rt2f,tfdes,rtm1f,rtm2f;
	float gammdeg,rtmf,vm,vm1,vm2,rtm1,rtm2,rtm,vtm1,vtm2,vc;
	float t,h,s,xnc,zemplos,zem1,zem2,rt1old,rt2old,vt1old,vt2old;
	float rm1old,rm2old,vm1old,vm2old,at1,at2,am1,am2;
	float atg,rt1k,rt2k,rm1k,rm2k,xncg,atplosg,rho,q,gamt;
	float tgo,atplos,xlam,xlamd;
	FILE *fptr;
	fptr=fopen("DATFIL.TXT","w");
	apn=0;
	xnp=3.;
    rt1=0.;
	rt2=200000.;
	rm1=170000.;
	rm2=0.;
	vt=6000.;
	rt2des=50000.;
	gamtdeg=45.;
	beta=500.;
	betest=500.;
	xnclimg=7.;
	xnclim=xnclimg*32.2;
	vt1=vt*cos(gamtdeg/57.3);
	vt2=-vt*sin(gamtdeg/57.3);
	initial(rt2des,rt1,rt2,vt1,vt2,betest,&rt1f,&rt2f,&tfdes);
	rtm1f=rt1f-rm1;
	rtm2f=rt2f-rm2;
	gammdeg=57.3*atan2(rtm2f,rtm1f);
	rtmf=sqrt(rtm1f*rtm1f+rtm2f*rtm2f);
	vm=rtmf/tfdes;
	vm1=vm*cos(gammdeg/57.3);
	vm2=vm*sin(gammdeg/57.3);
	rtm1=rt1-rm1;
	rtm2=rt2-rm2;
	rtm=sqrt(rtm1*rtm1+rtm2*rtm2);
	vtm1=vt1-vm1;
	vtm2=vt2-vm2;
	vc=-(rtm1*vtm1+rtm2*vtm2)/rtm;
 	t=0.;
	h=.01;
	s=0.;
	xnc=0.;
	zemplos=0.;
	zem1=0.;
	zem2=0.;
L5:
	if(vc<0.)
		goto L999;
 	if(rtm<1000.)
		h=.0002;
	else
		h=.01;
	s=s+h;
	rt1old=rt1;
	rt2old=rt2;
	vt1old=vt1;
	vt2old=vt2;
	rm1old=rm1;
	rm2old=rm2;
	vm1old=vm1;
	vm2old=vm2;
 	step=1;
	goto L200;
L66:	
	step=2;
	rt1=rt1+h*vt1;
	rt2=rt2+h*vt2;
	vt1=vt1+h*at1;
	vt2=vt2+h*at2;
	rm1=rm1+h*vm1;
	rm2=rm2+h*vm2;
	vm1=vm1+h*am1;
	vm2=vm2+h*am2;
 	t=t+h;
	goto L200;
L55:
	rt1=.5*(rt1old+rt1+h*vt1);
	rt2=.5*(rt2old+rt2+h*vt2);
	vt1=.5*(vt1old+vt1+h*at1);
	vt2=.5*(vt2old+vt2+h*at2);
	rm1=.5*(rm1old+rm1+h*vm1);
	rm2=.5*(rm2old+rm2+h*vm2);
	vm1=.5*(vm1old+vm1+h*am1);
	vm2=.5*(vm2old+vm2+h*am2);
 	if(s>=.09999){
		s=0.;
		atg=sqrt(at1*at1+at2*at2)/32.2;
		rt1k=rt1/1000.;
		rt2k=rt2/1000.;
		rm1k=rm1/1000.;
		rm2k=rm2/1000.;
		xncg=xnc/32.2;
		atplosg=atplos/32.2;
		vm=sqrt(vm1*vm1+vm2*vm2);
		printf("%10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n",
		t,rt1k,rt2k,rm1k,rm2k,atg,xncg,atplosg);
		fprintf(fptr,"%10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n",
		t,rt1k,rt2k,rm1k,rm2k,atg,xncg,atplosg);
 	}
	goto L5;
L200:
 	if(rt2<=30000.)
		rho=.002378*exp(-rt2/30000.);
	else
		rho=.0034*exp(-rt2/22000.);
	vt=sqrt(vt1*vt1+vt2*vt2);
	q=.5*rho*vt*vt;
	gamt=atan2(-vt2,vt1);
	at1=-32.2*q*cos(gamt)/beta;
	at2=-32.2+32.2*q*sin(gamt)/beta;
	rtm1=rt1-rm1;
	rtm2=rt2-rm2;
	rtm=sqrt(rtm1*rtm1+rtm2*rtm2);
	vtm1=vt1-vm1;;
	vtm2=vt2-vm2;;
	vc=-(rtm1*vtm1+rtm2*vtm2)/rtm;
	xlam=atan2(rtm2,rtm1);
	xlamd=(rtm1*vtm2-rtm2*vtm1)/(rtm*rtm);
	atplos=-at1*sin(xlam)+at2*cos(xlam);
	if(apn==1)
		xnc=xnp*vc*xlamd+.5*xnp*atplos;
	else
		xnc=xnp*vc*xlamd;
	if(xnc>xnclim)
		xnc=xnclim;
	if(xnc<-xnclim)
		xnc=-xnclim;
	am1=-xnc*sin(xlam);
	am2=xnc*cos(xlam);
     	if(step<=1)
		goto L66;
	else
		goto L55;
L999:
 	atg=sqrt(at1*at1+at2*at2)/32.2;
	rt1k=rt1/1000.;
	rt2k=rt2/1000.;
	rm1k=rm1/1000.;
	rm2k=rm2/1000.;
	xncg=xnc/32.2;
	atplosg=atplos/32.2;
	vm=sqrt(vm1*vm1+vm2*vm2);
	printf("%10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n",
		t,rt1k,rt2k,rm1k,rm2k,atg,xncg,atplosg);
	fprintf(fptr,"%10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n",
		t,rt1k,rt2k,rm1k,rm2k,atg,xncg,atplosg);
 	printf("%10.4f\n",rtm);
	fclose (fptr);
	return 0;
}
	
	void initial(float rt2des,float rt1,float rt2,float vt1,float vt2,float betest,
		float *rt1f,float *rt2f,float *tfdes)
{
	int step;
	float t,h,rt1old,rt2old,vt1old,vt2old,at1,at2,rho,vt,q;
	float gamt,beta;
	t=0.;
	h=.01;
	beta=betest;
L5:
	if(rt2<=rt2des)
		goto L999;
	rt1old=rt1;
	rt2old=rt2;
	vt1old=vt1;
	vt2old=vt2;
 	step=1;
	goto L200;
L66:
	step=2;
	rt1=rt1+h*vt1;
	rt2=rt2+h*vt2;
	vt1=vt1+h*at1;
	vt2=vt2+h*at2;
 	t=t+h;
	goto L200;
L55:
	rt1=.5*(rt1old+rt1+h*vt1);
	rt2=.5*(rt2old+rt2+h*vt2);
	vt1=.5*(vt1old+vt1+h*at1);
	vt2=.5*(vt2old+vt2+h*at2);
	goto L5;
L200:
 	if(rt2<=30000.)
		rho=.002378*exp(-rt2/30000.);
	else
		rho=.0034*exp(-rt2/22000.);
	vt=sqrt(vt1*vt1+vt2*vt2);
	q=.5*rho*vt*vt;
	gamt=atan2(-vt2,vt1);
	at1=-32.2*q*cos(gamt)/beta;
	at2=-32.2+32.2*q*sin(gamt)/beta;
     	if(step<=1)
		goto L66;
	else
		goto L55;
L999:
 	*rt1f=rt1;
	*rt2f=rt2;
	*tfdes=t;
}
	
	
