clc; clear all; close all; tic;

%This code used PN eq to get 3 types of evasive maneuvers : 
% First case [poly] :  we make the escape maneuver as nth degree polynomial
             % and we try to find the optimized solution using PN eq 
% Second case [Trap-general] :  we make the escape maneuver as trapezoidal
             % acceleration, and we sacn all the domain for t1 and t2 which lead to
             % otimized solution
% Third case [Tooth - square] : special case of trap - symmetric
%% Configuration Parameters

polyORtrapORtrapSymm = 'trap'; % choose which type of target acceleration
N=70; %number of runs
G = 32.2; %gravity

%=====Parameters of the polynomial function=====
degree = 3; %Polynomial degree

%=====Parameters of the trapezoid function=====
amax=12*G; %max accelartion
t0=0;  %beginning of the ramp
%t1  %end of the ramp & beginning of the const. acc.
%t2  %end of the const. acc. & beginning of the decent
t3=10;  %end of the decent

%---cost function weighting---
% f=w1*XNC + w2*RTM;
w1=1;
w2=10^5;

%% Prepareing the variables for each case : 

if strcmp(polyORtrapORtrapSymm,'poly')
    cs=zeros(N,degree); % Matrix of Coefficients: generating vector of zeros
    f=zeros(N,1); % cost array: generating vector of zeros
    
elseif strcmp(polyORtrapORtrapSymm,'trap')
    cs=zeros(N,N,5); % Matrix of coefficients trap ; 5 is length(c)
    f=zeros(N,N); % cost matrix: generating matrix of zeros
    
    t1 = zeros(N,1);
    t2 = zeros(N,N);
elseif strcmp(polyORtrapORtrapSymm,'trapSymm')
    cs=zeros(N,5); % Matrix of coefficients trap ; 5 is length(c)
    f=zeros(N,1); % cost array: generating vector of zeros
    
    t1 = zeros(N,1);
    t2 = zeros(N,1);
end

%% Three cases of evasive maneuvers 

if strcmp(polyORtrapORtrapSymm,'poly') %===== [1] POLY SOLUTION =====
    
    for(i=1:1:N)
        %making the escape manuver for the target as a nth degree polynomial
        c = 3* G* rand(1,degree); %random generation for the coeff. of the polynomial
        
        cs(i,:)=c;
        [RTM,XNCmax] = PNeq (polyORtrapORtrapSymm, degree, c ,0);
        
        %---cost function---
        % filling the vector with the values of the cost fn
        f(i,1)=w1*(XNCmax^2) + w2*(RTM^2)
        
    end
    
    [s,R]= max(f) %get index and value of max. f
    cs(R,:); %the numbers that generate max cost
    [RTM,XNCmax] = PNeq (polyORtrapORtrapSymm, degree, cs(R,:),1) %draw the optimized solution
    
 %-------------------------------------------------------------------------
    
elseif strcmp(polyORtrapORtrapSymm,'trap')  %=== [2] TRAPEZOIDAL SOLUTION ===
    %         c= rand(1,2); %generating array of two rundom number
    %         ct=sort(c); %rearange the array, starting from the smallest value
    %         t1= t3*ct(1); %t1>t2 & scaling with t3
    %         t2 = t3*ct(2);
    
    for(i=1:1:N+1)
        t1(i)= t0 + i*((t3-t0)/N);
        
        for (j=i:1:N+1)
            t2(i,j)= t1(i,1)+ j*((t3-t0)/N);
            c = [t0 t1(i,1) t2(i,j) t3 amax];
            cs(i,j,:)= c;
            [RTM,XNCmax] = PNeq (polyORtrapORtrapSymm, degree, c ,0);
            
            %---cost function---
            % filling the vector with the values of the cost fn
            f(i,j)= w1*(XNCmax^2) + w2*(RTM^2);
        end
        
    end
    
    %f(:); %put all the matrix elements in column array order
    [M,I] = max(f(:)); %get the max of this array and its index
    [I_row, I_col] = ind2sub(size(f),I); %refer the index to the original matrix
    
    %[s,R]= max(f)
    
    cs(I_row, I_col,:); %the numbers that generate max cost
    [RTM,XNCmax] = PNeq (polyORtrapORtrapSymm, degree, cs(I_row, I_col,:),1)
    %I_row
    %I_col
    
 %-------------------------------------------------------------------
    
elseif strcmp(polyORtrapORtrapSymm,'trapSymm')  %=== [3] Tooth - Square SOLUTION (speacial case of trap)===
    
    for(i=1:1:N+1/2)
        t1(i)= t0 + (i-1)*((t3-t0)/N);
        t2(i)= t3 - (i-1)*((t3-t0)/N);
        
        c = [t0 t1(i) t2(i) t3 amax];
        cs(i,:)= c;
        [RTM,XNCmax] = PNeq (polyORtrapORtrapSymm, degree, c ,0);
        
        %---cost function---
        % filling the vector with the values of the cost fn
        f(i,1)= w1*(XNCmax^2) + w2*(RTM^2);
        
        
    end
    
    [s,R]= max(f) %get index and value of max. f
    cs(R,:); %the numbers that generate max cost
    [RTM,XNCmax] = PNeq (polyORtrapORtrapSymm, degree, cs(R,:),1) %draw the optimized solution
end


%output=[ArrayT',ArrayRT1',ArrayRT2',ArrayRM1',ArrayRM2',ArrayXNCG',ArrayRTM' ];
%save datfil.txt output /ascii
disp '*** Simulation Complete'

toc;