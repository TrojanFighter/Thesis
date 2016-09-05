clc; clear all; close all; tic;

%TWO-DIMENSIONAL TACTICAL MISSILE-TARGET ENGAGEMENT SIMULATION
%solving using the second-order Runge–Kutta numerical integration technique

%% Configuration Parameters

polyORtrap = 'trap'; % choose which type of target acceleration
N=50; %number of runs
G = 32.2; %gravity


%=====Parameters of the polynomial function=====
degree = 2; %Polynomial degree

%=====Parameters of the trapezoid function=====
amax=12*G; %max accelartion
t0=0;  %beginning of the ramp
%t1 will be randumly generated %end of the ramp & beginning of the const. acc.
%t2 will be randumly generated %end of the const. acc. & beginning of the decent
t3=4;  %end of the decent


%---cost function weighting---
% f=w1*XNC + w2*RTM;
w1=1;
w2=10^5;


f=zeros(1,N); % cost array: generating vector of zeros

if strcmp(polyORtrap,'poly')
    cs=zeros(N,degree); % Matrix of Coefficients: generating vector of zeros
elseif strcmp(polyORtrap,'trap')
    cs=zeros(N,5); % Matrix of coefficients trap ; 5 is length(c) 
end

for(i=1:1:N)
    
    if strcmp(polyORtrap,'poly')
        %making the escape manuver for the target as a nth degree polynomial
        c = 3* G* rand(1,degree); %random generation for the coeff. of the polynomial
        
        
    elseif strcmp(polyORtrap,'trap')
        c= rand(1,2); %generating array of two rundom number
        ct=sort(c); %rearange the array, starting from the smallest value
        t1= t3*ct(1); %t1>t2 & scaling with t3 
        t2 = t3*ct(2);
        c = [t0 t1 t2 t3 amax]; 
    end
    
    cs(i,:)=c; 
    [RTM,XNCmax] = final (polyORtrap, degree, c ,0);
    
    %---cost function---
    % filling the vector with the values of the cost fn
    f(1,i)=w1*(XNCmax^2) + w2*(RTM^2)
    
    
end

[s,R]=max(f)
%R=find(max(f))
cs(R,:); %the numbers that generate max cost
[RTM,XNCmax] = final (polyORtrap, degree, cs(R,:),1)


%output=[ArrayT',ArrayRT1',ArrayRT2',ArrayRM1',ArrayRM2',ArrayXNCG',ArrayRTM' ];
%save datfil.txt output /ascii
disp '*** Simulation Complete'

toc;