#include <stdio.h>
#include <stdlib.h>
#include <math.h>
main()
{
	int step;
	float rt1,rt2,vt,gamtdeg,beta,vt1,vt2,t,h,s,rt1old,rt2old;
	float vt1old,vt2old,at1,at2,atg,rt1k,rt2k,rho;
	float q,gamt;
	FILE *fptr;
	fptr=fopen("DATFIL.TXT","w");
     	rt1=0.;
	rt2=100000.;
	vt=6000.;
	gamtdeg=45.;
	beta=500.;
	vt1=vt*cos(gamtdeg/57.3);
	vt2=-vt*sin(gamtdeg/57.3);
 	t=0.;
	h=.01;
	s=0.;
L5:
	if(rt2<0.)goto L999;
	s=s+h;
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
 	if(s>=.09999){
		s=0.;
		atg=sqrt(at1*at1+at2*at2)/32.2;
		rt1k=rt1/1000.;
		rt2k=rt2/1000.;
		vt=sqrt(vt1*vt1+vt2*vt2);
		atg=sqrt(at1*at1+at2*at2)/32.2;
		printf("%10.4f %10.4f %10.4f %10.4f %10.4f\n",t,rt1k,rt2k,atg,vt);
		fprintf(fptr,"%10.4f %10.4f %10.4f %10.4f %10.4f\n",t,rt1k,rt2k,atg,vt);
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
     	if(step<=1)
		goto L66;
	else
		goto L55;
L999:
 	atg=sqrt(at1*at1+at2*at2)/32.2;
	rt1k=rt1/1000.;
	rt2k=rt2/1000.;
	vt=sqrt(vt1*vt1+vt2*vt2);
	atg=sqrt(at1*at1+at2*at2)/32.2;
	printf("%10.4f %10.4f %10.4f %10.4f %10.4f\n",t,rt1k,rt2k,atg,vt);
	fprintf(fptr,"%10.4f %10.4f %10.4f %10.4f %10.4f\n",t,rt1k,rt2k,atg,vt);
	fclose (fptr);
	return 0;
}
