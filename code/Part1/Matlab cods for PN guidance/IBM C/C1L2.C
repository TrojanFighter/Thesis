#include <stdio.h>
#include <stdlib.h>
#include <math.h>
main()
{
	int n;
	float g,x,ts,y,t,ytheory;
	FILE *fptr;
	fptr=fopen("DATFIL.TXT","w");
	g=.5;
	x=1.;
	ts=.1;
	y=0.;
	t=0.;
	n=0;
	ytheory=1.-pow(1.-g,n);
	printf("%10.3f %10.3f %10.3f \n",t,y,ytheory);
	fprintf(fptr,"%10.3f %10.3f %10.3f \n",t,y,ytheory);
	for(n=1;n<=20;n=n+1){
		y=y+g*(x-y);
		t=n*ts;
		ytheory=1.-pow(1.-g,n);
		printf("%10.3f %10.3f %10.3f \n",t,y,ytheory);
		fprintf(fptr,"%10.3f %10.3f %10.3f \n",t,y,ytheory);
	}
	fclose (fptr);
	return 0;
}
