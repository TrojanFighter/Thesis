#include <stdio.h>
#include <stdlib.h>
#include <math.h>
main()
{
	int step;
	float h,vm,beta,t,s,gamdeg,vm1,vm2,rm1,rm2,rm1old,rm2old,vm1old,vm2old;
	float am1,am2,rm1k,rm2k,rho,q,gam,drag;
	FILE *fptr;
	fptr=fopen("DATFIL.TXT","w");
	h=.01;
	vm=3000.;
	beta=1000.;
	t=0.;
	s=0.;
	gamdeg=45.;
	vm1=vm*cos(gamdeg/57.3);
	vm2=vm*sin(gamdeg/57.3);
	rm1=0.;
	rm2=0.;
L10:
	if(t>0.&&rm2<=0.)goto L999;
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
	rm1k=rm1/1000.;
	rm2k=rm2/1000.;
	printf("%6.3f %6.3f %6.3f\n", t,rm1k,rm2k);
	fprintf(fptr,"%6.3f %6.3f %6.3f\n", t,rm1k,rm2k);
	goto L10;
 L200:
 	if(rm2<30000.)
		rho=.002378*exp(-rm2/30000);
	else
		rho=.0034*exp(-rm2/22000);
	vm=sqrt(vm1*vm1+vm2*vm2);
	q=.5*rho*vm*vm;
	gam=atan2(vm2,vm1);
	drag=q*32.2/beta;
 	am1=-drag*cos(gam);
	am2=-32.2-drag*sin(gam);
	if (step<2)
		goto L66;
	else
		goto L55;
L999:
	vm=sqrt(vm1*vm1+vm2*vm2);
	rm1k=rm1/1000.;
	rm2k=rm2/1000.;
	printf("%6.3f %6.3f %6.3f\n", t,rm1k,rm2k);
	fprintf(fptr,"%6.3f %6.3f %6.3f\n", t,rm1k,rm2k);
	fclose (fptr);
	return 0;
}
