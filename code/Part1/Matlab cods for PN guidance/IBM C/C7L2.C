#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
void gauss(float signoise,float *lamnoise);
main()
{
	clock_t start, secs;
	float z[1000],vc,nt,yic,vm,hedeg,beta,np,signoise,ts,t,s;
	float tf,z1,y,yd,ydic,h,gfilter,hfilter,lamh,lamdh,nc;
	float yold,ydold,ydd,lamnoise,lam,res,tgo,rtm,mean,sigma;
	int run,i,step,quit;
	FILE *fptr;
begin:
	
	fptr=fopen("DATFIL.TXT","w");
	quit=0;
	vc=4000.;
	nt=96.6;
	printf("current values for closing velocity and target acceleration %f\t%f\n",vc,nt);
	printf ("enter closing velocity and target acceleration:");
	scanf ("%f%f",&vc,&nt);
	printf("you entered %f\t%f\n",vc,nt);
	start=clock();
	yic=hedeg=0.;
	vm=3000.;
	beta=.8;
	np=3.;
	signoise=.001;
	ts=.1;
	run=50;
	for (tf=.5; tf<=10.; tf=tf+.5){
		z1=0.;
		for (i=1; i<=run; i=i+1){
			y=yic;
			yd=-vm*hedeg/57.3;
			ydic=yd;
			t=s=0.;
			h=.01;
			gfilter=1.-beta*beta;
			hfilter=(1.-beta)*(1.-beta);
			lamh=lamdh=nc=0.;
L10:
			if(t>(tf-.0001))
				goto L999;
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
			y=.5*(y+yold+h*yd);
			yd=.5*(yd+ydold+h*ydd);
			s=s+h;
			if(s<(ts-.0001))
				goto L10;
			s=0.;
			gauss(signoise,&lamnoise);
			res=lam-(lamh+ts*lamdh)+lamnoise;
			lamh=gfilter*res+lamh+ts*lamdh;
			lamdh=hfilter*res/ts+lamdh;
			nc=np*vc*lamdh;
			goto L10;
L200:		
			tgo=tf-t+.00001;
			rtm=vc*tgo;
			lam=y/(vc*tgo);
			ydd=nt-nc;
			if (step<2)
				goto L66;
			else
				goto L55;
L999:
			z[i]=y;
			z1=z[i]+z1;
			mean=z1/i;
		}
		sigma=0.;
		z1=0.;
		for (i=1; i<=run; i=i+1){
			z1=(z[i]-mean)*(z[i]-mean)+z1;
			if(i==1)
				sigma=0.;
			else
				sigma=sqrt(z1/(i-1));
		}
			
		printf("%6.3f %6.3f %6.3f\n", tf,sigma,mean);
		fprintf(fptr,"%6.3f %6.3f %6.3f\n", tf,sigma,mean);
	}
	fclose (fptr);
	secs=(clock()-start)/60.;
	printf( "That took %lu seconds.\n", secs );
	printf ("do you want to quit (0=no,1=yes):");
	scanf ("%d",&quit);
	printf("you entered %d\n",quit);
	if(quit==0)
	goto begin;
}
 	
	
 	void gauss(float signoise,float *lamnoise)
{
	float sum,x,y,temp;
	int j;
	sum=0.;
	for (j=1; j<=6; j=j+1){
		x=rand();
		y=(float)x/(float)RAND_MAX;
		sum=sum+y;
		}
	temp=sum-3.;
	*lamnoise=1.414*temp*signoise;
}

	
