#include <stdio.h>
#include <stdlib.h>
#include <math.h>
void gauss(float sig,float *x);
main()
{
	int step;
	float tau,phi,t,s,h,sig,y,x,yold,yd,sigplus,sigminus;
	FILE *fptr;
	fptr=fopen("DATFIL.TXT","w");
	tau=.2;
	phi=1.;
	t=0.;
	h=.01;
	sig=sqrt(phi/h);
	y=0.;
L10:	if(t>4.999)goto L999;
	gauss(sig,&x);
 	yold=y;
	step=1;
	goto L200;
L66:	step=2;
 	y=y+h*yd;
	t=t+h;
	goto L200;
L55:
 	y=(yold+y)/2.+.5*h*yd;
	sigplus=sqrt(phi*(1.-exp(-2.*t/tau))/(2.*tau));
	sigminus=-sigplus;
	printf("%6.3f %6.3f %6.3f %6.3f\n", t,y,sigplus,sigminus);
	fprintf(fptr,"%6.3f %6.3f %6.3f %6.3f\n", t,y,sigplus,sigminus);
	goto L10;
L200:
 	yd=(x-y)/tau;
	if (step<2)
		goto L66;
	else
		goto L55;
L999:
	fclose (fptr);
	return 0;
}

	void gauss(float sig,float *x)
{
	float sum,y,temp;
	int j,x1;
	sum=0.;
	for (j=1; j<=12; j=j+1){
		x1=rand();
		y=(float)x1/(float)RAND_MAX;
		sum=sum+y;
		}
	temp=sum-6.;
	*x=temp*sig;
}
