#include <stdio.h>
#include <stdlib.h>
#include <math.h>
main()
{
	int qswitch,step;
	float xnp,displace,thom,vm,vt,rm1,rm2,rt1,rt2,vt1,vt2,t,s;
	float rtm1,rtm2,rtm,xlam,vm1,vm2,vtm1,vtm2,vc,tgo,h;
	float rt1old,rt2old,rm1old,rm2old,vm1old,vm2old;
	float am1,am2,xnc,xlamd;
	FILE *fptr;
	fptr=fopen("DATFIL.TXT","w");
	xnp=3.;
	displace=200.;
	thom=1.;
	vm=3000.;
	vt=1000.;
	rm1=0.;
	rm2=1000.;
	rt1=20000.;
	rt2=1000.;
	qswitch=0;
	vt1=-vt;
	vt2=0.;
	t=0.;
	s=0.;
	rtm1=rt1-rm1;
	rtm2=rt2-rm2;
	rtm=sqrt(rtm1*rtm1+rtm2*rtm2);
	xlam=atan2(rtm2,rtm1);
	vm1=vm;
	vm2=0.;
	vtm1=vt1-vm1;
	vtm2=vt2-vm2;
	vc=-(rtm1*vtm1+rtm2*vtm2)/rtm;
	tgo=rtm/vc;
	h=.005;
 L10:
 	if(vc<0.)goto L999;
	tgo=rtm/vc;
	if(tgo<.3)h=.00005;
	if(tgo<=thom&&qswitch==0){
		qswitch=1;
		rt2=rt2+displace;
	}
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
 	rt1=rt1+h*vt1;
	rt2=rt2+h*vt2;
	rm1=rm1+h*vm1;
	rm2=rm2+h*vm2;
	vm1=vm1+h*am1;
	vm2=vm2+h*am2;
	t=t+h;
	goto L200;
L55:
 	rt1=.5*(rt1old+rt1+h*vt1);
	rt2=.5*(rt2old+rt2+h*vt2);
	rm1=.5*(rm1old+rm1+h*vm1);
	rm2=.5*(rm2old+rm2+h*vm2);
	vm1=.5*(vm1old+vm1+h*am1);
	vm2=.5*(vm2old+vm2+h*am2);
	s=s+h;
	if(s<.049999)goto L10;
	s=0.;
	printf("%10.3f %10.3f %10.3f %10.3f %10.3f %10.3f\n",t,rm1,rm2,rt1,rt2,xnc/32.2);
	fprintf(fptr,"%10.3f %10.3f %10.3f %10.3f %10.3f %10.3f\n",t,rm1,rm2,rt1,rt2,xnc/32.2);
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
	if(step<=1)
		goto L66;
	else
		goto L55;
L999:
 	printf("%10.3f %10.3f\n",t,rtm);
	fclose (fptr);
	return 0;
}
