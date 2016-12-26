#include <stdio.h>
#include <stdlib.h>
#include <math.h>
void pzadd(float a,float *b,float*c);
main()
{
	double f[5][5],x[5][5],fx[5][5],fxt[5][5],fxfxt[5][5];
	double xold[5][5],k0[5][5],q[5][5];
	double k1[5][5],k2[5][5],k3[5][5],xd[5][5];
	double a[2][5],at[5][2],ax[2][5],axat[2][2];
	int order,step,i,j,k;
	double t,tnew,s,h,np,tau,nt,vc,tf,tgo,phis,sigy,signl;
	
	FILE *fptr;
	fptr=fopen("DATFIL.TXT","w");
	order=4;
	t=s=0.;
	tnew=t;
	h=.01;
	np=3.;
	tau=1.;
	nt=96.6;
	vc=4000.;
	tf=10.;
	tgo=tf-t+.00001;
	phis=nt*nt/tf;
	for (i=1; i<=order; i=i+1){
		for (j=1; j<=order; j=j+1){
			f[i][j]=0.;
			x[i][j]=0.;
			q[i][j]=0.;
		}
	}
	f[1][2]=1.;
 	f[2][1]=-np/(tau*tgo);
	f[2][3]=1.;
	f[2][4]=np*vc/tau;
	f[4][1]=1./(tau*vc*tgo);
	f[4][4]=-1./tau;
	q[3][3]=phis;
L5:
	if (t>=tf)
		goto L999;
	for (i=1; i<=order; i=i+1){
		for (j=1; j<=order; j=j+1){
			xold[i][j]=x[i][j];
		}
	}
L20:
	step=1;
	goto L200;
L40:
	step=2;
	for (i=1; i<=order; i=i+1){
		for (j=1; j<=order; j=j+1){
			k0[i][j]=xd[i][j];
		}
	}
	tnew=t+.5*h;
	for (i=1; i<=order; i=i+1){
		for (j=1; j<=order; j=j+1){
			x[i][j]=xold[i][j]+.5*h*k0[i][j];
		}
	}
	goto L200;
L41:
	step=3;
	for (i=1; i<=order; i=i+1){
		for (j=1; j<=order; j=j+1){
			k1[i][j]=xd[i][j];
		}
	}
	tnew=t+.5*h;
	for (i=1; i<=order; i=i+1){
		for (j=1; j<=order; j=j+1){
			x[i][j]=xold[i][j]+.5*h*k1[i][j];
		}
	}
	goto L200;
L42:
	step=4;
	for (i=1; i<=order; i=i+1){
		for (j=1; j<=order; j=j+1){
			k2[i][j]=xd[i][j];
		}
	}
	tnew=t+h;
	for (i=1; i<=order; i=i+1){
		for (j=1; j<=order; j=j+1){
			x[i][j]=xold[i][j]+h*k2[i][j];
		}
	}
	goto L200;
L43:
	for (i=1; i<=order; i=i+1){
		for (j=1; j<=order; j=j+1){
			k3[i][j]=xd[i][j];
		}
	}
	t=tnew;
	for (i=1; i<=order; i=i+1){
		for (j=1; j<=order; j=j+1){
			x[i][j]=xold[i][j]+h*(k0[i][j]+2.*(k1[i][j]+k2[i][j])+k3[i][j])/6.;
		}
	}
	s=s+h;
	if (s<=.09999)
		goto L5;
	s=0.;
	a[1][1]=np/(tau*tgo);
	a[1][2]=0.;
	a[1][3]=0.;
	a[1][4]=-np*vc/tau;
	for (i=1; i<=1; i=i+1){
		for (j=1; j<=order; j=j+1){
			ax[i][j]=0.;
			for (k=1; k<=order; k=k+1){
				ax[i][j]=ax[i][j]+a[i][k]*x[k][j];
			}
		}
	}
	for (i=1; i<=1; i=i+1){
		for (j=1; j<=order; j=j+1){
			at[j][i]=a[i][j];
		}
	}
 	for (i=1; i<=1; i=i+1){
		for (j=1; j<=order; j=j+1){
			axat[i][j]=0.;
			for (k=1; k<=order; k=k+1){
				axat[i][j]=axat[i][j]+ax[i][k]*at[k][j];
			}
		}
	}
	sigy=sqrt(x[1][1]);
	signl=sqrt(axat[1][1]);
	printf("%6.3f %6.3f %6.3f\n", t,sigy,signl);
	fprintf(fptr,"%6.3f %6.3f %6.3f\n", t,sigy,signl);
	goto L5;
L200:
	tgo=tf-tnew+.00001;
	f[2][1]=-np/(tau*tgo);
	f[4][1]=1./(tau*vc*tgo);
	for (i=1; i<=order; i=i+1){
		for (j=1; j<=order; j=j+1){
			fx[i][j]=0.;
			for (k=1; k<=order; k=k+1){
				fx[i][j]=fx[i][j]+f[i][k]*x[k][j];
			}
		}
	}
	for (i=1; i<=order; i=i+1){
		for (j=1; j<=order; j=j+1){
			fxt[j][i]=fx[i][j];
		}
	}
	for (i=1; i<=order; i=i+1){
		for (j=1; j<=order; j=j+1){
			fxfxt[i][j]=fx[i][j]+fxt[i][j];
		}
	}
	for (i=1; i<=order; i=i+1){
		for (j=1; j<=order; j=j+1){
			xd[i][j]=fxfxt[i][j]+q[i][j];
		}
	}
	if(step==1)
		goto L40;
	if(step==2)
		goto L41;
	if(step==3)
		goto L42;
	else
		goto L43;
L999:
	fclose (fptr);
	return 0;
}
 	
