#include <stdio.h>
#include <stdlib.h>
#include <math.h>
void predict(double t,double rm1,double rm2,double vm1,double vm2,
		double *x,double xnc,double rt1,double rt2,double betah);
main()
{
	double h,vm,beta,betah,t,s,gamdeg,vm1,vm2,rm1,rm2,xnc,rt1,rt2;
	double rtm1,rtm2,rtm,xlim,rm1old,rm2old,vm1old,vm2old,x,x1;
	double delx,xncp,x2,dxdnc,delxnc,rm1k,rm2k,rt1k,rt2k,xncg,xlam;
	double xne1,xne2,rho,q,gam,drag,am1,am2,xlimg;
	int step;
	FILE *fptr1;
	fptr1=fopen("DATFIL.TXT","w");
	h=.01;
	vm=3000.;
	beta=1000.;
	betah=1000.;
	xlimg=5.;
	t=0.;
	s=0.;
	gamdeg=45.;
	vm1=vm*cos(gamdeg/57.3);
	vm2=vm*sin(gamdeg/57.3);
	rm1=0.;
	rm2=0.;
	xnc=0.;
	rt1=60000.;
	rt2=0.;
	rtm1=rt1-rm1;
	rtm2=rt2-rm2;
	rtm=sqrt(rtm1*rtm1+rtm2*rtm2);
	xlim=xlimg*32.2;
L10:
	if(t>0.&rm2<=0.)goto L999;
	if(rtm<1000.)
		h=.0001;
	else
		h=.01;
 	rm1old=rm1;
	rm2old=rm2;
	vm1old=vm1;
	vm2old=vm2;
	step=1;
	goto L200;
L66:
	step=2;
 	rm1=rm1+h*vm1;
	rm2=rm2+h*vm2;
	vm1=vm1+h*am1;
	vm2=vm2+h*am2;
	t=t+h;
	goto L200;
L55:
 	rm1=.5*(rm1old+rm1+h*vm1);
	rm2=.5*(rm2old+rm2+h*vm2);
 	vm1=.5*(vm1old+vm1+h*am1);
	vm2=.5*(vm2old+vm2+h*am2);
	s=s+h;
	if(s<.99999)goto L10;
	s=0.;
	if(t>30.){
		predict(t,rm1,rm2,vm1,vm2,&x,xnc,rt1,rt2,betah);
		x1=x;
		delx=rt1-x1;
		xncp=xnc+1.;
		predict(t,rm1,rm2,vm1,vm2,&x,xncp,rt1,rt2,betah);
		x2=x;
		dxdnc=(x2-x1)/(xncp-xnc);
		delxnc=delx/dxdnc;
		xnc=xnc+delxnc;
		if(xnc>xlim)xnc=xlim;
		if(xnc<-xlim)xnc=-xlim;
	}
	rm1k=rm1/1000.;
	rm2k=rm2/1000.;
	rt1k=rt1/1000.;
	rt2k=rt2/1000.;
	xncg=xnc/32.2;
	printf("%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n"
		,t,rm1k,rm2k,rt1k,rt2k,xncg);
	fprintf(fptr1,"%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n"
		,t,rm1k,rm2k,rt1k,rt2k,xncg);
	goto L10;
L200:
	rtm1=rt1-rm1;
	rtm2=rt2-rm2;
	rtm=sqrt(rtm1*rtm1+rtm2*rtm2);
	xlam=atan2(rtm2,rtm1);
	xne1=-xnc*sin(xlam);
	xne2=xnc*cos(xlam);
 	if(rm2<30000.)
		rho=.002378*exp(-rm2/30000);
	else
		rho=.0034*exp(-rm2/22000);
	vm=sqrt(vm1*vm1+vm2*vm2);
	q=.5*rho*vm*vm;
	gam=atan2(vm2,vm1);
	drag=q*32.2/beta;
 	am1=-drag*cos(gam)+xne1;
	am2=-32.2-drag*sin(gam)+xne2;
	if (step<2)
		goto L66;
	else
		goto L55;
L999:
	vm=sqrt(vm1*vm1+vm2*vm2);
	rm1k=rm1/1000.;
	rm2k=rm2/1000.;
	rt1k=rt1/1000.;
	rt2k=rt2/1000.;
	xncg=xnc/32.2;
	printf("%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n"
		,t,rm1k,rm2k,rt1k,rt2k,xncg);
	fprintf(fptr1,"%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n"
		,t,rm1k,rm2k,rt1k,rt2k,xncg);
	printf("%6.3f\n",rtm);
	fclose (fptr1);
	return 0;
}
	
void predict(double tp,double rm1p,double rm2p,double vm1p,
	double vm2p,double *rm1f,double xnc1f,double rt1p,double rt2p,double betah)
{
	int step;
	double h,rm1,rm2,vm1,vm2,xnc,rt1,rt2,beta,t,rtm1,rtm2,rtm;
	double rm1old,rm2old,vm1old,vm2old,xlam,xne1,xne2,rho;
	double vm,q,gam,drag,am1,am2;
	h=.01;
	rm1=rm1p;
	rm2=rm2p;
	vm1=vm1p;
	vm2=vm2p;
	xnc=xnc1f;
	rt1=rt1p;
	rt2=rt2p;
	beta=betah;
	t=tp;
	rtm1=rt1-rm1;
	rtm2=rt2-rm2;
	rtm=sqrt(rtm1*rtm1+rtm2*rtm2);
L10:
	if(rm2<0.)goto L999;
	if(rtm<1000.)
		h=.0001;
	else
		h=.01;
 	rm1old=rm1;
	rm2old=rm2;
	vm1old=vm1;
	vm2old=vm2;
	step=1;
	goto L200;
L66:
	step=2;
 	rm1=rm1+h*vm1;
	rm2=rm2+h*vm2;
	vm1=vm1+h*am1;
	vm2=vm2+h*am2;
	t=t+h;
	goto L200;
L55:
 	rm1=.5*(rm1old+rm1+h*vm1);
	rm2=.5*(rm2old+rm2+h*vm2);
 	vm1=.5*(vm1old+vm1+h*am1);
	vm2=.5*(vm2old+vm2+h*am2);
	goto L10;
L200:
	rtm1=rt1-rm1;
	rtm2=rt2-rm2;
	rtm=sqrt(rtm1*rtm1+rtm2*rtm2);
	xlam=atan2(rtm2,rtm1);
	xne1=-xnc*sin(xlam);
	xne2=xnc*cos(xlam);
 	if(rm2<30000.)
		rho=.002378*exp(-rm2/30000);
	else
		rho=.0034*exp(-rm2/22000);
	vm=sqrt(vm1*vm1+vm2*vm2);
	q=.5*rho*vm*vm;
	gam=atan2(vm2,vm1);
	drag=q*32.2/beta;
 	am1=-drag*cos(gam)+xne1;
	am2=-32.2-drag*sin(gam)+xne2;
	if (step<2)
		goto L66;
	else
		goto L55;
L999:
	*rm1f=rm1;
}
