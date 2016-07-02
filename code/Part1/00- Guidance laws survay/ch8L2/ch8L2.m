clc; clear all; close all; tic;

%adjoint simulation of the optimal single time constant guidance system (fig8.17-P.181)
%with various guidance law options
%solving using the second-order Runge–Kutta numerical integration technique

%% Simulation inputs

XNT=96.6;  % 3-g target accelration (target manuver)[error source]
XNP=4.;  %effictive navigation ratio
TAU=1.;  %guidance system time const
TF=10.;  %total flight time of the engagement
VM=3000.; %magnitude of the missile velocity
HEDEG=-20.; %heading error (in degrees)
APN=0;  %proportional navigation if APN = 0 
        %augmented proportional navigation if APN = 1
        %optimal guidance if APN = 2
T=0.;
S=0.;
TP=T+.00001;
%stats of the adjoint block diagram 
X1=0.;
X2=0.;
X3=1.; %I.C
X4=0.;
XNPP=0.;  %optimal effictive navigation ratio
H=.01;  %step size
HE=HEDEG/57.3;  %heading error (in rad)
n=0.;  %counter on points
%%
while TP<=(TF-1e-5) %program is stopped when the current time equals the flight time
	X1OLD=X1;
	X2OLD=X2;
	X3OLD=X3;
	X4OLD=X4;
   STEP=1;
   FLAG=0;
	while STEP<=1
      if FLAG==1
         STEP=2;
         X1=X1+H*X1D;
         X2=X2+H*X2D;
         X3=X3+H*X3D;
         X4=X4+H*X4D;
         TP=TP+H;
      end
      TGO=TP+.00001;
      if APN==0  %proportional navigation
         C1=XNP/(TGO*TGO);
         C2=XNP/TGO;
         C3=0.;
         C4=0.;
      elseif APN==1  %augmented proportional navigation
         C1=XNP/(TGO*TGO);
         C2=XNP/TGO;
         C3=.5*XNP;
         C4=0.;
      else %optimal guidance
         X=TGO/TAU;
         TOP=6.*X*X*(exp(-X)-1.+X);
         BOT1=2*X*X*X+3.+6.*X-6.*X*X;
         BOT2=-12.*X*exp(-X)-3.*exp(-2.*X);
         XNPP=TOP/(.0001+BOT1+BOT2); %optimal effictive navigation ratio
         C1=XNPP/(TGO*TGO);
         C2=XNPP/TGO;
         C3=.5*XNPP;
         C4=-XNPP*(exp(-X)+X-1.)/(X*X);
      end
      %adjoint differential equations (fig.8.17)
      X1D=X2+C3*X4/TAU;
      X2D=X3+C2*X4/TAU;
      X3D=C1*X4/TAU;
      X4D=-X4/TAU-X2+C4*X4/TAU;
      FLAG=1;
   end
   FLAG=0;
	X1=(X1OLD+X1)/2+.5*H*X1D;
	X2=(X2OLD+X2)/2+.5*H*X2D;
	X3=(X3OLD+X3)/2+.5*H*X3D;
    X4=(X4OLD+X4)/2+.5*H*X4D;
   S=S+H;
	if S>=.0999
      S=0.;
      n=n+1;
      ArrayTP(n)=TP;
      ArrayXMNT(n)=XNT*X1;
      ArrayXMHE(n)=-VM*HE*X2;
   end
end
figure
plot(ArrayTP,ArrayXMNT),grid
xlabel('Flight Time (Sec)')
ylabel('Target Maneuver Miss (Ft)')
clc
output=[ArrayTP',ArrayXMNT',ArrayXMHE'];
save datfil.txt output  -ascii
disp 'simulation finished'


toc;