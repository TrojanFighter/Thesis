clc; clear all; close all; tic;

%ADJOINT SIMULATION OF SINGLE-LAG GUIDANCE SYSTEM(fig8.11-P.175)
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
T=0.;  %time
S=0.;
TP=T+.00001;
%stats of the adjoint block diagram 
X1=0;
X2=0;
X3=1; %I.C
X4=0;
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
      
      %adjoint differential equations
      X1D=X2+X4*XNP*APN/(2.*TAU);  %miss distance sensivity due a step target manuver
      X2D=X3+XNP*X4/(TAU*TP);
      X3D=XNP*X4/(TAU*TP*TP);
      X4D=-X4/TAU-X2;
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
      ArrayXMNT(n)=XNT*X1;  %miss distance due a step target manuver
      ArrayXMHE(n)=-VM*HE*X2;  %miss distance due a step heading error
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