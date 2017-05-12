function output(solution,options,data)

%OUTPUT - Format and display the solution. Estimate the discretization error.
%
% Syntax:  output(solutions,options,data)
%
% Inputs:
%    solutions - Structure containing the solution
%    options - Display options
%    data - structure containing matrices to format data 
%
% Other m-files required: ppcreate.m
% Subfunctions: estimateError
% MAT-files required: none
%
% Copyright (C) 2010 Paola Falugi, Eric Kerrigan and Eugene van Wyk. All Rights Reserved.
% This code is published under the BSD License.
% Department of Electrical and Electronic Engineering,
% Imperial College London London  England, UK 
% ICLOCS (Imperial College London Optimal Control) 5 May 2010
% iclocs@imperial.ac.uk

%------------- BEGIN CODE --------------

% Define some variables
[nt,np,n,m,ng,nb,M,N,ns]=deal(data.sizes{:});
vdat=data.data;
t0=data.t0;

z=solution.z;
Vx=data.map.Vx;
Vu=data.map.Vu;
f=data.functions{3};

% Generate time vector
if nt; tf=z(1);else   tf=data.tf(1);end

T=(tf-t0)*[0;cumsum(data.tau)]*data.Nm/ns+data.k0;
Tdec=T(1:ns:end);

% Extract design parameters
if np; 
p=solution.p';
P=repmat(p,M,1);
else P=[];
p=[];
end

% Extract states and controls

X=reshape(Vx*z,n,M)';
usp=reshape(data.map.Vu*z,m,N)';
U=kron(usp,ones((M-1)/N,1));
U=[U;U(end,:)];


% Define piecewise polynomials
F=f(X,U,P,T,vdat);
U(end,:)=[];
Xp=cell(n,1);dXp=cell(n,1);Up=cell(m,1);

for i=1:n % Cubic Hermite interpolation
 [Xp{i}, dXp{i}]=Hsplines(T,X(:,i),F(:,i));
end

for i=1:m % Piecewise constant polynomials
    Up{i}=mkpp(T,U(:,i)');
end


% Display computation time
if (options.print.time)
    disp('computation time:');disp(solution.computation_time);
end

% Display minimized cost
if (options.print.cost)
    disp('minimized cost:');disp(costFunction(z,data));
end

% Plot states
if (options.plot.states)
    figure(1)
    for i=1:n
        hold on
        plot(T,ppval(Xp{i},T),'r')
        
    end
    title('Optimal state trajectories')
    xlabel('Time')
    grid on; axis tight;
end


% Plot inputs
if (options.plot.inputs)
    figure(2)
    for i=1:m
        hold on
        plot(T,ppval(Up{i},T),'r')
    end
    title('Optimal input sequences')
    xlabel('Time')
    grid on; axis tight;
end

% Plot multipliers
if (options.plot.multipliers==1);

  
  switch(data.options.NLPsolver)
   case{'ipopt'} 
     % Estimate the adjoint variables
       lambda_midpoint=reshape(solution.multipliers.lambda(n+1:n*M),n,M-1)';
   case{'fmincon'}  
      lambda_midpoint=reshape(solution.multipliers.eqnonlin(n+1:n*M),n,M-1)';
  end    
  lambda_midpoint=lambda_midpoint(ns:ns:M-1,:);
  if M<=2
      lambda=lambda_midpoint;
    else
      lambda=interp1((1.5:(M-1)/ns+1)',lambda_midpoint,(1:(M-1)/ns+1)','linear','extrap');
   end     
  
   
    figure(3)
    for i=1:n
        hold on
        plot(Tdec,lambda(:,i),'r')
    end
    title('Adjoint variables')
    xlabel('Time')
    grid on; axis tight;



end

if options.print.relative_local_error
  Error=estimateError(Xp,Up,p,dXp,tf,n,m,f,M,ns,data);
  figure(4)
    hold on
    plot(T(2:end),Error)
    xlabel('Time')
    title('Relative local error') 
    grid on; axis tight; 
end
%------------- END OF CODE --------------




function Error=estimateError(Xp,Up,p,dXp,tf,n,m,f,M,ns,data)

%ESTIMATEERROR - Estimate the discretization error
%Represent the solution by piecewise cubic polynomials p(t) using 
%Hermite interpolation and intergrate dp/dt-f(p) using Romberg quadrature
%
% Syntax:  Error=estimateError(Xp,Up,p,dXp,tf,n,m,f,M,ns,data)
%
% Inputs:
%    Xp  - Structure containing polynomial coefficients for states.
%    Up  - Structure containing polynomial coefficients for inputs.
%    dXp - Structure containing polynomial coefficients for state derivatives.
%
% Output:
%    Error - Global discretization error
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
%------------- BEGIN CODE --------------


q=11;% estimated relative error 10^(-q);
T=(tf-data.t0)*[0;cumsum(data.tau)]*data.Nm/ns+data.k0;
vdat=data.data;

%Computation of scale weights for the discretization error
wi=zeros(n,1);
for i=1:n
     wi(i)=max(max(abs(ppval(Xp{i},T)),abs(ppval(dXp{i},T))));
end

Error=zeros(M-1,1);         %Pre-allocation

for k=1:M-1


a=T(k); b=T(k+1);
h = 2.^((1:q)-1)*(b-a)/2^(q-1);             % These are the intervals used.
k1 = 2.^((q-2):-1:-1)*2+1;                  % Index into the intervals.
tq=a:h(1):b;

Xq=zeros(length(tq),n);
dXq=zeros(length(tq),n);
for i=1:n
    xp=ppval(Xp{i},tq);
    Xq(:,i)=xp(:);
    dxp=ppval(dXp{i},tq);
    dXq(:,i)=dxp(:); 
end

Uq=zeros(length(tq),m);
for i=1:m
    up=ppval(Up{i},tq); 
    Uq(:,i)=up(:);
end

% Extract design parameters if specified and convert to cells
if ~isempty(p) 
    P=repmat(p,k1(1),1);
else
  P=[];  
end

F=dXq-f(Xq,Uq,P,tq',vdat);           % Function evaluations.


% computation of the absolute local error on step k


for kk=1:n
  R = zeros(1,q);                                           % Pre-allocation.
  romberg=zeros(n,1);
  % Define the starting vector:
  for ii = 1:q
	R(ii) = 0.5*h(ii)*(F(1,kk)+2*...
                       sum(F(k1(end-ii+1):k1(end-ii+1)-1:end-1,kk))+F(end,kk));
  end
% Interpolations:
for jj = 2:q
    jpower = (4^(jj-1)-1);
    for ii = 1:(q-jj+1)
        R(ii) = R(ii)+(R(ii)-R(ii+1))/jpower; 
    end 
 end
romberg(kk) = R(1);
end

Error(k)=max(romberg./(wi+1));
end
%------------- END OF CODE --------------