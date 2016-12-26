#include <stdio.h>
#include <stdlib.h>
#include <math.h>
main()
{
	int step;
	float xisp1,xisp2,xmf1,xmf2,wpay,delv,delv1,delv2,amax1,amax2;
	float top2,bot2,wp2,ws2,wtot2,trst2,tb2,top1,bot1,wp1,ws1;
	float wtot,trst1,tb1,delvk,h,t,s,v,vold,a;
	float ag,vk,wgt,trst;
	FILE *fptr;
	fptr=fopen("DATFIL.TXT","w");
	xisp1=250.;
	xisp2=250.;
	xmf1=.85;
	xmf2=.85;
	wpay=100.;
	delv=20000.;
	delv1=.3333*delv;
	delv2=.6667*delv;
	amax1=10.;
	amax2=10.;
	top2=wpay*(exp(delv2/(xisp2*32.2))-1.);
	bot2=1/xmf2-((1.-xmf2)/xmf2)*exp(delv2/(xisp2*32.2));
	wp2=top2/bot2;
	ws2=wp2*(1-xmf2)/xmf2;
	wtot2=wp2+ws2+wpay;
	trst2=amax2*(wpay+ws2);
	tb2=xisp2*wp2/trst2;
	top1=wtot2*(exp(delv1/(xisp1*32.2))-1.);
	bot1=1/xmf1-((1.-xmf1)/xmf1)*exp(delv1/(xisp1*32.2));
	wp1=top1/bot1;
	ws1=wp1*(1-xmf1)/xmf1;
	wtot=wp1+ws1+wtot2;
	trst1=amax1*(wtot2+ws1);
	tb1=xisp1*wp1/trst1;
	delvk=delv/1000.;
	h=.01;
	t=0.;
	s=0.;
	v=0.;
L10:	if(t>(tb1+tb2))goto L999;
 	vold=v;
	step=1;
	goto L200;
L66:	step=2;
	v=v+h*a;
	t=t+h;
	goto L200;
L55:	v=(vold+v)/2+.5*h*a;
 	s=s+h;
	if(s<.99999)goto L10;
 	s=0.;
	ag=a/32.2;
	vk=v/1000.;
	printf("%8.2f %8.2f %8.2f\n",t,vk,ag);
	fprintf(fptr,"%8.2f %8.2f %8.2f\n",t,vk,ag);
	goto L10;
L200:
 	if(t<tb1){
		wgt=-wp1*t/tb1+wtot;
		trst=trst1;
	}
	else if(t<(tb1+tb2)){
		wgt=-wp2*t/tb2+wtot2+wp2*tb1/tb2;
		trst=trst2;
	}
	else{
		wgt=wpay;
		trst=0.;
	}
	a=32.2*trst/wgt;
	if(step<=1)
		goto L66;
	else
		goto L55;
L999:
 	ag=a/32.2;
	vk=v/1000.;
 	printf("%8.2f %8.2f %8.2f\n",t,vk,ag);
	fprintf(fptr,"%8.2f %8.2f %8.2f\n",t,vk,ag);
	fclose (fptr);
	return 0;
}
