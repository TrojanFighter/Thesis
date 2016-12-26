#include <stdio.h>
#include <stdlib.h>
#include <math.h>
void gauss(double signoise,double *xnoise);
main()
{
	double phi[4][4],p[4][4],m[4][4],phip[4][4],phipphit[4][4],gain[4][2];
	double q[4][4],hmat[2][4],hm[2][4],mht[4][2];
	double phit[4][4];
	double phthpr[4][4],hmht[2][2],ht[4][2],kh[4][4],idn[4][4],ikh[4][4];
	double hmhtr,hmhtrinv,xnoise,rt2db,vt2db,res,erry,sp11,errv;
 	double sp22,errbeta,sp33,rt2k;
	int order,step,i,j,k;
	double rt2,vt2,beta,rt2h,vt2h,betah,ts,tf,q33,t,s,h,signoise;
	double rt2old,vt2old,rt2d,vt2d,rhoh,f21,f22,f23;
	FILE *fptr1;
	FILE *fptr2;
	fptr1=fopen("DATFIL.TXT","w");
	fptr2=fopen("COVFIL.TXT","w");
	rt2=100000.;
	vt2=-6000.;
	beta=500.;
	rt2h=100025.;
	vt2h=-6150.;
	betah=800.;
	order=3;
	ts=.05;
	tf=30.;
	q33=0./tf;
	t=0.;
	s=0.;
	h=.001;
	signoise=sqrt(500.);
	for (i=1; i<=order; i=i+1){
		for (j=1; j<=order; j=j+1){
			phi[i][j]=0.;
			p[i][j]=0.;
			q[i][j]=0.;
			idn[i][j]=0.;
		}
	}
 	idn[1][1]=1.;
	idn[2][2]=1.;
	idn[3][3]=1.;
	p[1][1]=signoise*signoise;
	p[2][2]=20000.;
	p[3][3]=300.*300.;
	for (i=1; i<=order; i=i+1){
		hmat[1][i]=0.;
		ht[i][1]=0.;
	}
	hmat[1][1]=1.;
	ht[1][1]=1.;
L10:	if(rt2<0.)goto L999;
	rt2old=rt2;
	vt2old=vt2;
 	step=1;
	goto L200;
L66:	step=2;
 	rt2=rt2+h*rt2d;
	vt2=vt2+h*vt2d;
 	t=t+h;
	goto L200;
L55:
 	rt2=.5*(rt2old+rt2+h*rt2d);
	vt2=.5*(vt2old+vt2+h*vt2d);
 	s=s+h;
	if(s<=(ts-.00001))goto L10;
	s=0.;
	rhoh=.0034*exp(-rt2h/22000.);
	f21=-32.2*rhoh*vt2h*vt2h/(2.*22000.*betah);
	f22=rhoh*32.2*vt2h/betah;
	f23=-rhoh*32.2*vt2h*vt2h/(2.*betah*betah);
	phi[1][1]=1.;
	phi[1][2]=ts;
	phi[2][1]=f21*ts;
	phi[2][2]=1.+f22*ts;
	phi[2][3]=f23*ts;
	phi[3][3]=1.;
	q[2][2]=f23*f23*q33*ts*ts*ts/3.;
	q[2][3]=f23*q33*ts*ts/2.;
	q[3][2]=q[2][3];
	q[3][3]=q33*ts;
	for (i=1; i<=order; i=i+1){
		for (j=1; j<=order; j=j+1){
			phit[j][i]=phi[i][j];
		}
	}
	for (i=1; i<=order; i=i+1){
		for (j=1; j<=order; j=j+1){
			phip[i][j]=0.;
			for (k=1; k<=order; k=k+1){
				phip[i][j]=phip[i][j]+phi[i][k]*p[k][j];
			}
		}
	}
 	for (i=1; i<=order; i=i+1){
		for (j=1; j<=order; j=j+1){
			phipphit[i][j]=0.;
			for (k=1; k<=order; k=k+1){
				phipphit[i][j]=phipphit[i][j]+phip[i][k]*phit[k][j];
			}
		}
	}
	for (i=1; i<=order; i=i+1){
		for (j=1; j<=order; j=j+1){
			m[i][j]=phipphit[i][j]+q[i][j];
		}
	}
 	for (i=1; i<=1; i=i+1){
		for (j=1; j<=order; j=j+1){
			hm[i][j]=0.;
			for (k=1; k<=order; k=k+1){
				hm[i][j]=hm[i][j]+hmat[i][k]*m[k][j];
			}
		}
	}
	for (i=1; i<=1; i=i+1){
		for (j=1; j<=1; j=j+1){
			hmht[i][j]=0.;
			for (k=1; k<=order; k=k+1){
				hmht[i][j]=hmht[i][j]+hm[i][k]*ht[k][j];
			}
		}
	}
	hmhtr=hmht[1][1]+signoise*signoise;
	hmhtrinv=1./hmhtr;
	for (i=1; i<=order; i=i+1){
		for (j=1; j<=1; j=j+1){
			mht[i][j]=0.;
			for (k=1; k<=order; k=k+1){
				mht[i][j]=mht[i][j]+m[i][k]*ht[k][j];
			}
		}
	}
	for (i=1; i<=order; i=i+1){
		gain[i][1]=mht[i][1]*hmhtrinv;
	}
	for (i=1; i<=order; i=i+1){
		for (j=1; j<=order; j=j+1){
			kh[i][j]=0.;
			for (k=1; k<=1; k=k+1){
				kh[i][j]=kh[i][j]+gain[i][k]*hmat[k][j];
			}
		}
	}
	for (i=1; i<=order; i=i+1){
		for (j=1; j<=order; j=j+1){
			ikh[i][j]=idn[i][j]-kh[i][j];
		}
	}
	for (i=1; i<=order; i=i+1){
		for (j=1; j<=order; j=j+1){
			p[i][j]=0.;
			for (k=1; k<=order; k=k+1){
				p[i][j]=p[i][j]+ikh[i][k]*m[k][j];
			}
		}
	}
 	gauss(signoise,&xnoise);
	rt2db=vt2h;
	vt2db=.0034*32.2*vt2h*vt2h*exp(-rt2h/22000.)/(2.*betah)-32.2;
	res=rt2+xnoise-(rt2h+rt2db*ts);
	rt2h=rt2h+rt2db*ts+gain[1][1]*res;
	vt2h=vt2h+vt2db*ts+gain[2][1]*res;
	betah=betah+gain[3][1]*res;
	erry=rt2-rt2h;
	sp11=sqrt(p[1][1]);
	errv=vt2-vt2h;
	sp22=sqrt(p[2][2]);
	errbeta=beta-betah;
	sp33=sqrt(p[3][3]);
	rt2k=rt2/1000.;
	printf("%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n"
			,t,rt2k,rt2,rt2h,vt2,vt2h,beta,betah);
	fprintf(fptr1,"%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n"
			,t,rt2k,rt2,rt2h,vt2,vt2h,beta,betah);
	fprintf(fptr2,"%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n"
			,t,rt2k,erry,sp11,-sp11,errv,sp22,-sp22,errbeta,sp33,-sp33);
     	goto L10;
L200:
	rt2d=vt2;
	vt2d=.0034*32.2*vt2*vt2*exp(-rt2/22000.)/(2.*beta)-32.2;
	if (step<2)
		goto L66;
	else
		goto L55;
L999:
 	fclose (fptr1);
 	fclose (fptr2);
	return 0;
}

	void gauss(double signoise,double *xnoise)
{
	double sum,x,y,temp;
	int j;
	sum=0.;
	for (j=1; j<=6; j=j+1){
		x=rand();
		y=(float)x/(float)RAND_MAX;
		sum=sum+y;
		}
	temp=sum-3.;
	*xnoise=1.414*temp*signoise;
}
