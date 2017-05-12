function [problem,guess] = BatchFermentor
%BatchFermentor - Optimal control of a fed-batch ethanol fermentor
%
% Syntax:  [problem,guess] = BatchFermentor
%
% Outputs:
%    problem - Structure with information on the optimal control problem
%    guess   - Guess for state, control and multipliers.
%
% Other m-files required: none
% Subfunctions: stageCost, boundaryCost, f, g, b
% MAT-files required: none
%
% Copyright (C) 2010 Paola Falugi, Eric Kerrigan and Eugene van Wyk. All Rights Reserved.
% This code is published under the BSD License.
% Department of Electrical and Electronic Engineering,
% Imperial College London London  England, UK 
% ICLOCS (Imperial College London Optimal Control) 5 May 2010
% iclocs@imperial.ac.uk

%------------- BEGIN CODE --------------


% Initial time. t0<tf
problem.time.t0=0;


% Final time. Let tf_min=tf_max if tf is fixed.
problem.time.tf_min=126;     
problem.time.tf_max=126; 
guess.tf=126;

% Parameters bounds. pl=< p <=pu
problem.parameters.pl=[];
problem.parameters.pu=[];
guess.parameters=[];


% Initial conditions for system.
problem.states.x0=[1.5 0.0 0.0 7]; 


% Initial conditions for system. Bounds if x0 is free s.t. x0l=< x0 <=x0u
problem.states.x0l=[1.5 0.0 0.0 7]; 
problem.states.x0u=[1.5 0.0 0.0 7]; 

% State bounds. xl=< x <=xu
problem.states.xl=[0 0.00 0.00 0.0];
problem.states.xu=[40 50 25 10];

% Terminal state bounds. xfl=< xf <=xfu
problem.states.xfl=[0 0.00 0.00 0.00]; 
problem.states.xfu=[40 50 25 10];

% Guess the state trajectories with [x0 xf]
guess.states(:,1)=[1.5 30];
guess.states(:,2)=[0 8.5];
guess.states(:,3)=[0 0];
guess.states(:,4)=[7 10];


% Number of control actions (default: ZOH) 
% Set M= 0 if M = # of integration steps + 1
% Set M=-1 for a PWL approximation
problem.inputs.N=500;       
      
% Input bounds
problem.inputs.ul=[0];
problem.inputs.uu=[50];

% Guess the input sequences with [u0 uf]
guess.inputs(:,1)=[2, 10];        


% Choose the set-points if required
problem.setpoints.states=[];
problem.setpoints.inputs=[];

% Bounds for path constraint function gl =< g(x,u,p,t) =< gu
problem.constraints.gl=[];
problem.constraints.gu=[];

% Bounds for boundary constraints bl =< b(x0,xf,u0,uf,p,t0,tf) =< bu
problem.constraints.bl=[];
problem.constraints.bu=[];

% store the necessary problem parameters used in the functions
problem.data=[];

% Get function handles and return to Main.m
problem.functions={@L,@E,@f,@g,@b};

%------------- END OF CODE --------------

function stageCost=L(x,xr,u,ur,p,t,data)
% L - Returns the stage cost
% For direct transcription methods the function must be vectorized and
% xi, ui are column vectors taken as x(:,i) and u(:,i)
% 
% Syntax:  stageCost = L(x,xr,u,ur,p,t)
%
% Inputs:
%    x  - state vector
%    xr - state reference
%    u  - input
%    ur - input reference
%    p  - parameter
%    t  - time
%
% Output:
%    stageCost - Scalar stage cost
%
%------------- BEGIN CODE --------------


stageCost =0.00001*u(:,1).*u(:,1);
 

%------------- END OF CODE --------------


function boundaryCost=E(x0,xf,u0,uf,p,tf,data) 
% E - Returns the boundary value cost
%
% Syntax:  boundaryCost=E(x0,xf,u0,uf,p,tf)
%
% Inputs:
%    x0  - state at t=0
%    xf  - state at t=tf
%    u0  - input at t=0
%    uf  - input at t=tf
%    p   - parameter
%    tf  - final time
%
% Output:
%    boundaryCost - Scalar boundary cost
%
%------------- BEGIN CODE --------------

boundaryCost=-xf(2)*xf(4);

%------------- END OF CODE --------------


function dx = f(x,u,p,t,data)
% f - Returns the ODE right hand side where x'= f(x,u,p,t)
% For direct transcription methods: the function must be vectorized and
% xi, ui, pi are column vectors taken as x(:,i), u(:,i) and p(:,i). Each
% state corresponds to one column of dx.
% 
% 
% Syntax:  dx = f(x,u,p,t)
%
% Inputs:
%    x  - state vector
%    u  - input
%    p  - parameter
%    t  - time
%
% Output:
%    dx - time derivative of x
%
%------------- BEGIN CODE --------------

x1 = x(:,1);x2 = x(:,2);x3 = x(:,3);x4 = x(:,4);
u1 = u(:,1);

h1 = 0.11*(x3./(0.006*x1+x3));
h2 = 0.0055*(x3./(0.0001+x3.*(1+10*x3)));

dx(:,1) = (h1.*x1-u1.*(x1./500./x4));
dx(:,2) = (h2.*x1-0.01*x2-u1.*(x2./500./x4));
dx(:,3) = (-h1.*x1/0.47-h2.*x1/1.2-x1.*(0.029*x3./(0.0001+x3))+u1./x4.*(1-x3/500));
dx(:,4) = u1/500;

%------------- END OF CODE --------------


function c=g(x,u,p,t,data)
% g - Returns the path constraint function where gl =< g(x,u,p,t) =< gu
% For direct transcription methods the function must be vectorized and
% xi, ui, pi are column vectors taken as x(:,i), u(:,i) and p(:,i). Each
% constraint corresponds to one column of c
% 
% Syntax:  c=g(x,u,p,t)
%
% Inputs:
%    x  - state vector
%    u  - input
%    p  - parameter
%    t  - time
%
% Output:
%    c - constraint function
%
%------------- BEGIN CODE --------------

c=[];

%------------- END OF CODE --------------

function bc=b(x0,xf,u0,uf,p,tf,data)
% b - Returns the boundary constraints: bl =< bf(x0,xf,u0,uf,p,t0,tf) =< bu
%
% Syntax:  bc=b(x0,xf,u0,uf,p,tf)
%
% Inputs:
%    x0  - state at t=0
%    xf  - state at t=tf
%    u0  - input at t=0
%    uf  - input at t=tf
%    p   - parameter
%    tf  - final time
%
% Output:
%    bc -  evaluation of boundary function 
%
%------------- BEGIN CODE --------------

bc=[];
%------------- END OF CODE --------------






