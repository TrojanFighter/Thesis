#include <stdio.h>
#include <stdlib.h>
#include <math.h>
main()
{
	int step;
	float vm,vt,xnt,rm1ic,rm2ic,rt1ic,rt2ic,hedeg,xnp,ts;
	float rt1,rt2,rm1,rm2,betat,vt1,vt2,he,xnc,t,s;
	float rtm1,rtm2,rtm,thett,vm1,vm2,vtm1,vtm2,vc;
	float betatold,rt1old,rt2old,rm1old,rm2old,vm1old,vm2old; 
	float betatd,am1,am2,h,rt1k,rt2k,rm1k,rm2k;
	float thetm,rt,rm,xlam,xlamd;
	FILE *fptr;
	fptr=fopen("DATFIL.TXT","w");
	vm=3000.;
	vt=1000.;
	xnt=0.;
	rm1ic=0.;
	rm2ic=1.;
	rt1ic=40000.;
	rt2ic=10000.;
	hedeg=0.;
	xnp=10.;
	ts=.1;
	rt1=rt1ic;
	rt2=rt2ic;
	rm1=rm1ic;
	rm2=rm2ic;
	betat=0.;
	vt1=-vt*cos(betat);
	vt2=vt*sin(betat);
	he=hedeg/57.3;
	xnc=0.;
	t=0.;
	s=0.;
	rtm1=rt1-rm1;
	rtm2=rt2-rm2;
	rtm=sqrt(rtm1*rtm1+rtm2*rtm2);
	thett=atan2(rt2,rt1);
	vm1=vm*cos(thett+he);
	vm2=vm*sin(thett+he);
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
	rt1k=rt1/1000.;
	rt2k=rt2/1000.;
	rm1k=rm1/1000.;
	rm2k=rm2/1000.;
	printf("%8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n",t,rt1,rt2,rm1,rm2,xnc,rtm);
	fprintf(fptr,"%8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n",t,rt1,rt2,rm1,rm2,xnc,rtm);
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
	xnc=xnp*rm*(thett-thetm);
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
 	printf("%8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n",t,rt1,rt2,rm1,rm2,xnc,rtm);
	fclose (fptr);
	return 0;
}
