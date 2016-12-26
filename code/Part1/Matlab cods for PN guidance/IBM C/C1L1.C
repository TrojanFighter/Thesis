#include <stdio.h>
#include <stdlib.h>
#include <math.h>
main()
{
	int step;
	float w,t,s,y,yd,x,h,yold,ydold,ydd,ytheory;
	FILE *fptr;
	fptr=fopen("DATFIL.TXT","w");
	w=20.;
	t=0.;
	s=0.;
	y=0.;
	yd=0.;
	x=1.;
	h=.001;
L10:
	if(t>1.)goto L999;
	s=s+h;
	yold=y;
	ydold=yd;
	step=1;
	goto L200;
L66:
	step=2;
	y=y+h*yd;
	yd=yd+h*ydd;
	t=t+h;
	goto L200;
L55:
	y=.5*(yold+y+h*yd);
	yd=.5*(ydold+yd+h*ydd);
	if(s<.009999)goto L10;
	s=0.;
	ytheory=(1.-cos(w*t))/w;
	printf("%10.3f %10.3f %10.3f \n",t,y,ytheory);
	fprintf(fptr,"%10.3f %10.3f %10.3f \n",t,y,ytheory);
	goto L10;
L200:
	ydd=w*x-w*w*y;
	if(step<=1)
		goto L66;
	else
		goto L55;
L999:
	fclose (fptr);
	return 0;
}
