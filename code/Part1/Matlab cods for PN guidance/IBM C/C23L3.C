#include <stdio.h>
#include <stdlib.h>
#include <math.h>
main()
{
	int i;
	float zact,wact,k3,ta,zaf,waf,kr,w,xmag1,xmag2,xmag3;
	float gain,phase1,phase2,phase3,phase;
	FILE *fptr;
	fptr=fopen("DATFIL.TXT","w");
	zact=.7;
	wact=150.;
	k3=-1.89;
	ta=.457;
	zaf=.058;
	waf=25.3;
	kr=.1;
	for(i=2;i<=160;i=i+1){
		w=pow(10,.025*i-1);
		xmag1=sqrt(1+(w*ta)*(w*ta));
		xmag2=sqrt((1-(w/waf)*(w/waf))*(1-(w/waf)*(w/waf))+(2*zaf*w/waf)*(2*zaf*w/waf));
		xmag3=sqrt((1-(w/wact)*(w/wact))*(1-(w/wact)*(w/wact))+(2*zact*w/wact)*(2*zact*w/wact));
		gain=20*log10(-k3*kr*xmag1/(xmag2*xmag3));
		phase1=57.3*atan2(w*ta,1.);
		phase2=57.3*atan2(2*zaf*w/waf,1-(w/waf)*(w/waf));
		phase3=57.3*atan2(2*zact*w/wact,1-(w/wact)*(w/wact));
		phase=phase1-phase2-phase3;
		printf("%10.3f %10.3f %10.3f \n",w,gain,phase);
		fprintf(fptr,"%10.3f %10.3f %10.3f \n",w,gain,phase);
	}
 	fclose (fptr);
	return 0;
}
