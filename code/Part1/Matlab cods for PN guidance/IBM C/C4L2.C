#include <stdio.h>
#include <stdlib.h>
#include <math.h>
void gauss(float sig,float *xnoise);
main()
{
	int i,x,n,j,bin,k;
	float y,sum,xmax,xmin,range,tmp,pdf,ab,th,sig,xnoise;
	float h[2000],z[2000];
	FILE *fptr;
	sig=1.;
	n=1000;
	xmax=6.;
	xmin=-6.;
	bin=50;
	range=xmax-xmin;
	tmp=1./sqrt(6.28);
	fptr=fopen("DATFIL.TXT","w");
	for (i=1; i<=n; i=i+1){
		gauss(sig,&xnoise);
		z[i]=xnoise;
	}
	for (i=1; i<=bin; i=i+1){
		h[i]=0.;
	}
	for (i=1; i<=n; i=i+1){
		k=(((z[i]-xmin)/range)*bin)+.99;
		if(k<1)
			k=1;
		if(k>bin)
			k=bin;
		h[k]=h[k]+1;
	}
	for (k=1; k<=bin; k=k+1){
		pdf=(h[k]/n)*bin/range;
		ab=xmin+k*range/bin;
		th=tmp*exp(-ab*ab/2.);
		printf("%6.3f %6.3f %6.3f\n", ab,pdf,th);
		fprintf(fptr,"%6.3f %6.3f %6.3f\n", ab,pdf,th);
	}
	fclose (fptr);
	return 0;
}

	void gauss(float sig,float *xnoise)
{
	float sum,x,y,temp;
	int j;
	sum=0.;
	for (j=1; j<=6; j=j+1){
		x=rand();
		y=(float)x/(float)RAND_MAX;
		sum=sum+y;
		}
	temp=sum-3.;
	*xnoise=1.414*temp*sig;
}
	
 	
 	
 	
