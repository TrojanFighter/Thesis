clc; clear all; close all; tic;
%LINEARIZED ENGAGEMENT SIMULATION
%solving using the second-order Runge–Kutta numerical integration technique

%% Simulation inputs

VC=4000.; %closing velocity
XNT=0.; %target accelration (target manuver)[first error source]
Y=0.;
VM=3000.; %magnitude of the missile velocity
HEDEG=-20.; %heading error (in degrees)     [second error source]
TF=10.; %total flight time of the engagement
XNP=3.5; %effictive navigation ratio
YD=-VM*HEDEG/57.3;
T=0.; %time
H=.01; %step size (kept fixed for the entire flight)
S=0.;
n=0.; %counter on points
%%
while T<=(TF-1e-5) %program is stopped when the current time equals the flight time
   YOLD=Y;
   YDOLD=YD;   
   STEP=1;
   FLAG=0;
   while STEP<=1
      if FLAG==1
         STEP=2;
         Y=Y+H*YD;
         YD=YD+H*YDD;
         T=T+H;
      end
      TGO=TF-T+.00001; %time to go until the end of the flight
      XLAMD=(Y+YD*TGO)/(VC*TGO*TGO); %missile lead angle
      XNC=XNP*VC*XLAMD; %acceleration command(PN eq)
      YDD=XNT-XNC;
      FLAG=1;
   end
   FLAG=0;
   Y=.5*(YOLD+Y+H*YD);
   YD=.5*(YDOLD+YD+H*YDD);
   S=S+H;
   if S>=.0999
      S=0.;
      n=n+1;
      ArrayT(n)=T;
      ArrayY(n)=Y;
      ArrayYD(n)=YD;
      ArrayXNCG(n)=XNC/32.2;
   end
end
figure
plot(ArrayT,ArrayXNCG),grid
xlabel('Time (Sec)')
ylabel('Missile Acceleration (G)')
clc
output=[ArrayT',ArrayY',ArrayYD',ArrayXNCG'];
save datfil.txt output  -ascii
disp 'simulation finished'

toc;
