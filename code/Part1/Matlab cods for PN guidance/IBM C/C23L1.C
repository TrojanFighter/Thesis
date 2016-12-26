#include <stdio.h>
#include <stdlib.h>
#include <math.h>
main()
{
	int step;
	float vm,xncg,del,alt,a,diam,fr,xl,ctw,crw,hw,ctt,crt;
	float ht,xn,xcg,xhl,wgt,rho,swing,stail,sref,xlp,splan;
	float xcpn,an,ab,xcpb,xcpw,xmach,xiyy,tmp1,tmp2;
	float tmp3,tmp4,b,q,xnl;
	float p1,y1,y2,y3,y4,y5,y6,p2,p3,alftr,deltr,cna,cnd,cmap;
	float cma,cmd,xma,xmd,za,zd,wz,waf,zaf,xk1,xk2,ta,xk3,xkdc;
	float e,ed,t,h,s,eold,edold,edd;
	FILE *fptr;
	fptr=fopen("DATFIL.TXT","w");
	vm=3000.;
	xncg=10.;
	alt=0.;
	a=1000.;
 	diam=1.;
	fr=3.;
	xl=20.;
	ctw=0.;
	crw=6.;
	hw=2.;
	ctt=0.;
	crt=2.;
	ht=2.;
	xn=4.;
	xcg=10.;
	xhl=19.5;
	wgt=1000.;
	if(alt<=30000.)
		rho=.002378*exp(-alt/30000.);
	else
		rho=.0034*exp(-alt/22000.);
	swing=.5*hw*(ctw+crw);
	stail=.5*ht*(ctt+crt);
	sref=3.1416*diam*diam/4.;
	xlp=fr*diam;
	splan=(xl-xlp)*diam+1.33*xlp*diam/2.;
	xcpn=2*xlp/3;
	an=.67*xlp*diam;
	ab=(xl-xlp)*diam;
	xcpb=(.67*an*xlp+ab*(xlp+.5*(xl-xlp)))/(an+ab);
	xcpw=xlp+xn+.7*crw-.2*ctw;
	xmach=vm/a;
	xiyy=wgt*(3*((diam/2)*(diam/2))+xl*xl)/(12*32.2);
	tmp1=(xcg-xcpw)/diam;
	tmp2=(xcg-xhl)/diam;
	tmp3=(xcg-xcpb)/diam;
	tmp4=(xcg-xcpn)/diam;
	b=sqrt(xmach*xmach-1);
	q=.5*rho*vm*vm;
	p1=wgt*xncg/(q*sref);
	y1=2.+8*swing/(b*sref)+8*stail/(b*sref);
	y2=1.5*splan/sref;
	y3=8*stail/(b*sref);
	y4=2*tmp4+8*swing*tmp1/(b*sref)+8*stail*tmp2/(b*sref);
	y5=1.5*splan*tmp3/sref;
	y6=8*stail*tmp2/(b*sref);
	p2=y2-y3*y5/y6;
	p3=y1-y3*y4/y6;
	alftr=(-p3+sqrt(p3*p3+4.*p2*p1))/(2.*p2);
	deltr=-y4*alftr/y6-y5*alftr*alftr/y6;
	cna=2+1.5*splan*alftr/sref+8*swing/(b*sref)+8*stail/(b*sref);
	cnd=8*stail/(b*sref);
	cmap=2*tmp4+1.5*splan*alftr*tmp3/sref+8*swing*tmp1/(b*sref);
	cma=cmap+8*stail*tmp2/(b*sref);
	cmd=8*stail*tmp2/(b*sref);
	xma=q*sref*diam*cma/xiyy;
	xmd=q*sref*diam*cmd/xiyy;
	za=-32.2*q*sref*cna/(wgt*vm);
	zd=-32.2*q*sref*cnd/(wgt*vm);
	wz=sqrt((xma*zd-za*xmd)/zd);
	waf=sqrt(-xma);
	zaf=.5*waf*za/xma;
	xk1=-vm*(xma*zd-xmd*za)/(1845*xma);
	xk2=xk1;
	ta=xmd/(xma*zd-xmd*za);
	xk3=1845*xk1/vm;
	xkdc=1./xk1;
	e=0.;
	ed=0.;
	t=0;
	h=.0001;
	s=0;
L10:
	if (t>1.99999) goto L999;
	s=s+h;
	eold=e;
	edold=ed;
	step=1;
	goto L200;
L66:
	step=2;
	e=e+h*ed;
	ed=ed+h*edd;
	t=t+h;
	goto L200;
L55:
	e=.5*(eold+e+h*ed);
	ed=.5*(edold+ed+h*edd);
	if(s<=.0099999)goto L10;
	s=0.;
	printf("%10.3f %10.3f %10.3f \n",t,xnl,xncg);
	fprintf(fptr,"%10.3f %10.3f %10.3f \n",t,xnl,xncg);
	goto L10;
L200:
 	del=xkdc*xncg;
 	edd=waf*waf*(del-e-2.*zaf*ed/waf);
 	xnl=xk1*(e-edd/(wz*wz));
	if(step<=1)
		goto L66;
	else
		goto L55;
L999:
	fclose (fptr);
	return 0;
}
