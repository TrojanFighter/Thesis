#include <stdio.h>
#include <stdlib.h>
#include <math.h>
main()
{
	int qswitch,step;
	float xnp,tau,displace,thom,vm,vt,rm1,rm2,rt1,rt2;
	float vt1,vt2,t,s,rtm1,rtm2,xlam,vm1,vm2,vtm1,vtm2;
	float vc,tgo,xlamh,h,am1,am2,xlamhd,xnc;
	float rt1old,rt2old,rm1old,rm2old,vm1old,vm2old,xlamhold;
	float x,theory,rtm,rtmp;
	FILE *fptr;
	fptr=fopen("DATFIL.TXT","w");
	xnp=3.;
	tau=1.;
	displace=200.;
	for(thom=.1;thom<=5.;thom=thom+.1){
		vm=3000.;
		vt=1000.;
		rm1=0.;
		rm2=0.;
		rt1=20000.;
		rt2=0.;
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
		xlamh=0.;
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
		xlamhold=xlamh;
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
		xlamh=xlamh+h*xlamhd;
		t=t+h;
		goto L200;
L55:
 		rt1=.5*(rt1old+rt1+h*vt1);
		rt2=.5*(rt2old+rt2+h*vt2);
		rm1=.5*(rm1old+rm1+h*vm1);
		rm2=.5*(rm2old+rm2+h*vm2);
		vm1=.5*(vm1old+vm1+h*am1);
		vm2=.5*(vm2old+vm2+h*am2);
		xlamh=.5*(xlamhold+xlamh+h*xlamhd);
		s=s+h;
		goto L10;
L200:
 		rtm1=rt1-rm1;
		rtm2=rt2-rm2;
		rtm=sqrt(rtm1*rtm1+rtm2*rtm2);
		vtm1=vt1-vm1;
		vtm2=vt2-vm2;
		vc=-(rtm1*vtm1+rtm2*vtm2)/rtm;
		xlam=atan2(rtm2,rtm1);
		xlamhd=(xlam-xlamh)/tau;
		xnc=xnp*vc*xlamhd;
		am1=-xnc*sin(xlam);
		am2=xnc*cos(xlam);
		if(step<=1)
			goto L66;
		else
			goto L55;
L999:
 		x=thom/tau;
 		theory=displace*exp(-x)*(1.-2.*x+.5*x*x);
		if(rtm2>0.)
			rtmp=rtm;
		else
			rtmp=-rtm;
		printf("%10.3f %10.3f %10.3f \n",thom,rtmp,theory);
		fprintf(fptr,"%10.3f %10.3f %10.3f \n",thom,rtmp,theory);
 	}
 	fclose (fptr);
	return 0;
}
