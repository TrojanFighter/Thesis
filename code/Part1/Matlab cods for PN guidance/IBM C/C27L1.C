#include <stdio.h>
#include <stdlib.h>
#include <math.h>
main()
{	
	int step,apn;
	double ts,xlimg,vm,beta,t,s,gamdeg,vm1,vm2,rm1,rm2,xnc,rt1,rt2,rtm1,rtm2,rtm;
	double vtm1,vtm2,vc,xlim,rm1old,rm2old,vm1olfd,vm2old,rm1k,rm2k;
	double rt1k,rt2k,drag1,drag2,dragplos,atplos;
	double xncg,dragplosg,xlam,xlamd,rho,q,gam,drag,xne1,xne2,am1,am2,h,vm1old;
	FILE *fptr1;
	fptr1=fopen("DATFIL.TXT","w");
	ts=.1;
	apn=0;
	xlimg=5.;
	vm=3000.;
	beta=1000.;
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
	vtm1=0.-vm1;
	vtm2=0.-vm2;
	vc=-(rtm1*vtm1+rtm2*vtm2)/rtm;
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
	if(s<(ts-.00001))goto L10;
	s=0.;
	rm1k=rm1/1000.;
	rm2k=rm2/1000.;
	rt1k=rt1/1000.;
	rt2k=rt2/1000.;
	drag1=-drag*cos(gam);
	drag2=-drag*sin(gam)-32.2;
	dragplos=-drag1*sin(xlam)+drag2*cos(xlam);
	atplos=0.;
	if(t>30.){
		if(apn==0){
			xnc=3.*vc*xlamd;
		}
		else{
			xnc=3.*vc*xlamd+1.5*(atplos-dragplos);
		}
	}
	else{
		xnc=0.;
	}
	if(xnc>xlim)xnc=xlim;
	if(xnc<-xlim)xnc=-xlim;
	xncg=xnc/32.2;
	dragplosg=dragplos/32.2;
	printf("%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n"
		,t,rm1k,rm2k,rt1k,rt2k,xncg);
	fprintf(fptr1,"%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n"
		,t,rm1k,rm2k,rt1k,rt2k,xncg);
	goto L10;
L200:
	rtm1=rt1-rm1;
	rtm2=rt2-rm2;
	rtm=sqrt(rtm1*rtm1+rtm2*rtm2);
	vtm1=0.-vm1;
	vtm2=0.-vm2;
	vc=-(rtm1*vtm1+rtm2*vtm2)/rtm;
	xlam=atan2(rtm2,rtm1);
	xlamd=(rtm1*vtm2-rtm2*vtm1)/(rtm*rtm);
 	if(rm2<30000.){
		rho=.002378*exp(-rm2/30000);
	}
	else{
		rho=.0034*exp(-rm2/22000);
	}
	vm=sqrt(vm1*vm1+vm2*vm2);
	q=.5*rho*vm*vm;
	gam=atan2(vm2,vm1);
	drag=q*32.2/beta;
	xne1=-xnc*sin(xlam);
	xne2=xnc*cos(xlam);
 	am1=-drag*cos(gam)+xne1;
	am2=-32.2-drag*sin(gam)+xne2;
	if (step<2)
		goto L66;
	else
		goto L55;
L999:
	rm1k=rm1/1000.;
	rm2k=rm2/1000.;
	rt1k=rt1/1000.;
	rt2k=rt2/1000.;
	xncg=xnc/32.2;
	dragplosg=dragplos/32.2;
	printf("%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n"
		,t,rm1k,rm2k,rt1k,rt2k,xncg);
	fprintf(fptr1,"%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n"
		,t,rm1k,rm2k,rt1k,rt2k,xncg);
	printf("%6.3f\n",rtm);
	fclose (fptr1);
 	return 0;
}
	
	
