#include <stdio.h>
#include <stdlib.h>
#include <math.h>
main()
{
	int step;
	float vm,del,alt,a,diam,fr,xl,ctw,crw,hw,ctt,crt;
	float ht,xn,xcg,xhl,wgt,rho,swing,stail,sref,xlp,splan;
	float xcpn,an,ab,xcpb,xcpw,xmach,xiyy,tmp1,tmp2;
	float tmp3,tmp4,b,q,thd,alf,t,h,s,xnl;
	float thdold,alfold,thdd,alfd,xnlg,alfdeg,cn,cm;
	FILE *fptr;
	fptr=fopen("DATFIL.TXT","w");
	vm=3000.;
	del=5./57.3;
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
	thd=0;
	alf=0;
	t=0;
	h=.0025;
	s=0.;
L10:
	if(t>1.99999) goto L999;
	s=s+h;
	thdold=thd;
	alfold=alf;
	step=1;
	goto L200;
L66:
	step=2;
	thd=thd+h*thdd;
	alf=alf+h*alfd;
	t=t+h;
	goto L200;
L55:
	thd=.5*(thdold+thd+h*thdd);
	alf=.5*(alfold+alf+h*alfd);
	if(s<=.0099999)goto L10;
	s=0.;
	xnlg=xnl/32.2;
	alfdeg=alf*57.3;
	printf("%10.3f %10.3f %10.3f \n",t,xnlg,alfdeg);
	fprintf(fptr,"%10.3f %10.3f %10.3f \n",t,xnlg,alfdeg);
	goto L10;
L200:
 	cn=2*alf+1.5*splan*alf*alf/sref+8*swing*alf/(b*sref)+8*stail*(alf+del)/(b*sref);
	cm=2*alf*tmp4+1.5*splan*alf*alf*tmp3/sref+8*swing*alf*tmp1/(b*sref)+8*stail*(alf+del)*tmp2/(b*sref);
	thdd=q*sref*diam*cm/xiyy;
	xnl=32.2*q*sref*cn/wgt;
	alfd=thd-xnl/vm;
	if(step<=1)
		goto L66;
	else
		goto L55;
L999:
	fclose (fptr);
	return 0;
}
