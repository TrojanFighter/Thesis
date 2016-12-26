#include <stdio.h>
#include <stdlib.h>
#include <math.h>
main()
{
	int step;
	float fr,diam,xl,ctw,crw,hw,ctt,crt,ht,xn,xcg,xhl;
	float wact,zact,tf,vm,xncg,wcr,zeta,tau,alt,a,rho;
	float wgt,xnllin,xacc,swing,stail,sref,xlp,splan,xcpn;
	float an,ab,xcpb,xcpw,xmach,xiyy,tmp1,tmp2,tmp3,tmp4;
	float b,q,p1,y1,y2,y3,y4,y5,y6,p2,p3,alftr,deltr;
	float cna,cnd,cmap,cma,cmd,xma,xmd,za,zd,wz,waf;
	float zaf,xk1,xk2,ta,xk3,w,w0,z0,xkc,xka,xk0,xk,wi,xkr,xkdc;
	float edd,deldd,xd,xnl,thd,delc;
	float e,ed,deld,del,x,t,h,s,eold,edold,delold,deldold,xold;
	FILE *fptr;
	fptr=fopen("DATFIL","w");
	fr=3.;
	diam=1.;
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
	wact=150.;
	zact=.7;
	tf=1.;
	vm=3000.;
	xncg=10.;
	wcr=50.;
	zeta=.7;
	tau=.3;
	alt=0.;
	a=1000.;
	if(alt<=30000.)
		rho=.002378*exp(-alt/30000.);
	else
		rho=.0034*exp(-alt/22000.);
	wgt=1000.;
	xnllin=0.;
	xacc=xcg;
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
	y1=2+8*swing/(b*sref)+8*stail/(b*sref);
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
	w=(tau*wcr*(1+2.*zaf*waf/wcr)-1)/(2*zeta*tau);
	w0=w/sqrt(tau*wcr);
	z0=.5*w0*(2*zeta/w+tau-waf*waf/(w0*w0*wcr));
	xkc=(-w0*w0/(wz*wz)-1.+2.*z0*w0*ta)/(1.-2.*z0*w0*ta+w0*w0*ta*ta);
	xka=xk3/(xk1*xkc);
	xk0=-w*w/(tau*waf*waf);
	xk=xk0/(xk1*(1+xkc));
	wi=xkc*ta*w0*w0/(1+xkc+w0*w0/(wz*wz));
	xkr=xk/(xka*wi);
	xkdc=1.+1845./(xka*vm);
	e=0.;
	ed=0.;
	deld=0.;
	del=0.;
	x=0.;
	t=0;
	h=.0001;
	s=0;
L10:
	if(t>(tf-.00001))goto L999;
	s=s+h;
	eold=e;
	edold=ed;
	delold=del;
	deldold=deld;
	xold=x;
	step=1;
	goto L200;
L66:
	step=2;
	e=e+h*ed;
	ed=ed+h*edd;
	del=del+h*deld;
	deld=deld+h*deldd;
	x=x+h*xd;
	t=t+h;
	goto L200;
L55:
	e=.5*(eold+e+h*ed);
	ed=.5*(edold+ed+h*edd);
	del=.5*(delold+del+h*deld);
	deld=.5*(deldold+deld+h*deldd);
	x=.5*(xold+x+h*xd);
	if(s<=.0099999)goto L10;
	s=0.;
	printf("%10.3f %10.3f %10.3f \n",t,xnl,xncg);
	fprintf(fptr,"%10.3f %10.3f %10.3f \n",t,xnl,xncg);
	goto L10;
L200:
	thd=xk3*(e+ta*ed);
 	delc=xkr*(x+thd);
 	deldd=wact*wact*(delc-del-2.*zact*deld/wact);
 	edd=waf*waf*(del-e-2.*zaf*ed/waf);
 	xnl=xk1*(e-edd/(wz*wz));
 	xd=wi*(thd+xka*(xnl-xncg*xkdc));
	if(step<=1)
		goto L66;
	else
		goto L55;
L999:
	fclose (fptr);
	return 0;
}
