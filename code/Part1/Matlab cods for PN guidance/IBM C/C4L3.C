#include <stdio.h>
#include <stdlib.h>
#include <math.h>
main()
{
	int i,x,n,j;
	float y,sum,z[101],k,z1,x1,xmean,sigma;
	FILE *fptr;
	fptr=fopen("DATFIL.TXT","w");
	n=100;
	z1=0.;
	for (i=1; i<=n; i=i+1){
		sum=0.;
		for (j=1; j<=12; j=j+1){
			x=rand();
			y=(float)x/(float)RAND_MAX;
			sum=sum+y;
			x1=sum-6.;
		}
		z[i]=x1;
		z1=z[i]+z1;
		xmean=z1/i;
	}
	sigma=0.;
	z1=0.;
	for (i=1; i<=n; i=i+1){
		z1=(z[i]-xmean)*(z[i]-xmean)+z1;
		if(i==1)
			sigma=0.;
		else
			sigma=sqrt(z1/(i-1.));
		k=i;
		printf("%6.0f %6.3f\n", k,sigma);
		fprintf(fptr,"%6.0f %6.3f\n", k,sigma);
	}
	fclose (fptr);
	return 0;
}
 	
 	
 	
