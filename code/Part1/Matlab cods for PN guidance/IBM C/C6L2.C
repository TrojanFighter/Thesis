#include <stdio.h>
main()
{
	float vc=4000.;
	float nt=32.2;
	float yic=0.;
	float vm=3000.;
	float hedeg=0.;
	float tau=.5;
	float np=3.;
	float ta=10.;
	float r=-.01;
	int step;
	FILE *fptr;
	float y,yd,ydic,nl,elamdh,x4,x5,th,thh,t,h,s,tf;
	float yold,ydold,nlold,elamdhold,x4old,x5old,thold,thhold;
	float ydd,nld,elamdhd,x4d,x5d,thd,thhd,xnc;
	float tgo,lam,eps,dd,nc;
	fptr=fopen("DATFIL.TXT","w");
	for (tf=.1; tf<=10.; tf=tf+.1){;
		y=yic;
		yd=-vm*hedeg/57.3;
		ydic=yd;
		nl=0.;
		elamdh=0.;
		x4=0.;
		x5=0.;
		th=0.;
		thh=0.;
		t=0.;
		h=.01;
		s=0.;
ten:
		if (t>(tf-.0001))
			goto nine;
		yold=y;
		ydold=yd;
		nlold=nl;
		elamdhold=elamdh;
		x4old=x4;
		x5old=x5;
		thold=th;
		thhold=thh;
		step=1;
		goto two;
six:
		step=2;
		y=y+h*yd;
		yd=yd+h*ydd;
 		nl=nl+h*nld;
		elamdh=elamdh+h*elamdhd;
		x4=x4+h*x4d;
		x5=x5+h*x5d;
		th=th+h*thd;
		thh=thh+h*thhd;
		t=t+h;
		goto two;
five:
		y=.5*(yold+y+h*yd);
		yd=.5*(ydold+yd+h*ydd);
 		nl=.5*(nlold+nl+h*nld);
		elamdh=.5*(elamdhold+elamdh+h*elamdhd);
		x4=.5*(x4old+x4+h*x4d);
		x5=.5*(x5old+x5+h*x5d);
		th=.5*(thold+th+h*thd);
		thh=.5*(thhold+thh+h*thhd);
		s=s+h;
		goto ten;
		
two:
		tgo=tf-t+.00001;
		lam=y/(vc*tgo);
		eps=lam-th-thh+r*thh;
		dd=5.*eps/tau;
		elamdhd=5.*(dd-elamdh)/tau;
		nc=np*vc*elamdh;
		x4d=5.*(nc-x4)/tau;
		x5d=5.*(x4-x5)/tau;
		nld=5.*(x5-nl)/tau;
		thd=nl/vm+ta*nld/vm;
		thhd=dd-thd;
		ydd=nt-nl;
		if (step<2)
			goto six;
		else
			goto five;
nine:
		printf("%6.1f %6.2f \n",tf,y);
		fprintf(fptr,"%6.1f %6.2f \n",tf,y);
	}
		fclose (fptr);
		return 0;
}

