clc; clear all; close all; tic;
%Problem of finding the step response of a second order system 
%using Runge-Kutta

%%
W=20.; %natural frequancy of the second order system
T=0.;  %time
S=0.;

%initial condition
Y=0.;
YD=0.;

X=1.; %step input
H=.001; %step size(in time)
n=0.; %counter on points
%%

while T<=(1.-1e-5)  %looping over time duration
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
      YDD=W*X-W*W*Y; %second order D.eq of the system
      FLAG=1;
   end
   FLAG=0;
   
   %Runge-Kutta formula
   Y=.5*(YOLD+Y+H*YD);
   YD=.5*(YDOLD+YD+H*YDD);
   
   S=S+H;
   if S>=.000999
      S=0.;
      n=n+1;
      ArrayT(n)=T;
      ArrayY(n)=Y;
   end
end
figure
plot(ArrayT,ArrayY),grid
xlabel('Time (Sec)')
ylabel('y')
clc
output=[ArrayT',ArrayY'];
save datfil.txt output  -ascii
disp 'simulation finished'

toc;