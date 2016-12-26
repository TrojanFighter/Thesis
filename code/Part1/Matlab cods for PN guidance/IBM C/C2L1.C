#include <stdio.h>
#include <stdlib.h>
#include <math.h>
main()
{
	float vm,vt,xnt,hedeg,xnp,rm1,rm2,rt1,rt2;
	float beta,vt1,vt2,he,t,s,rtm1,rtm2,rtm,xlam,xlead;
	float thet,vm1,vm2,vtm1,vtm2,vc;
	float betaold,rt1old,rt2old,rm1old,rm2old,vm1old,vm2old;
	float betad,am1,am2,xlamd,xnc,h;
	int step;
	FILE *fptr;
	fptr=fopen("DATFIL.TXT","w");
	vm=3000.;
	vt=1000.;
	xnt=0.;
	hedeg=-20.;
	xnp=4.;
	rm1=0.;
	rm2=10000.;
	rt1=40000.;
	rt2=10000.;
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
		h=.0002;
	else
		h=.01;
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
	printf("%10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n",t,rt1,rt2,rm1,rm2,xnc);
	fprintf(fptr,"%10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n",t,rt1,rt2,rm1,rm2,xnc);
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
	xnc=xnp*vc*xlamd;
	am1=-xnc*sin(xlam);
	am2=xnc*cos(xlam);
	vt1=-vt*cos(beta);
	vt2=vt*sin(beta);
	betad=xnt/vt;
	if(step<=1)
		goto L66;
	else
		goto L55;
L999:
 	printf("%10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n",t,rt1,rt2,rm1,rm2,xnc);
	printf("%10.4f\n",rtm);
	fclose (fptr);
	return 0;
}
