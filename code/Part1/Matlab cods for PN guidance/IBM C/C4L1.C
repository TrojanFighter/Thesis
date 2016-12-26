#include <stdio.h>
#include <stdlib.h>
#include <math.h>
main()
{
	int i,x,n,j;
	float y,sum,x1,k;
	FILE *fptr;
	fptr=fopen("DATFIL.TXT","w");
	n=100;
	for (i=1; i<=n; i=i+1){
		sum=0.;
		for (j=1; j<=12; j=j+1){
			x=rand();
			y=(float)x/(float)RAND_MAX;
			sum=sum+y;
		}
		x1=sum-6.;
		k=i;
		printf("%6.0f %6.3f\n",k,x1);
		fprintf(fptr,"%6.0f %6.3f\n",k,x1);
	}
	fclose (fptr);
	return 0;
}
