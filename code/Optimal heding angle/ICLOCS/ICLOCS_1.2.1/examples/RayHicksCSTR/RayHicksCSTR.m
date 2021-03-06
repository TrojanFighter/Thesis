function [problem,guess] = RayHicksCSTR

%RayHicksCSTR- Define the optimal control problem for  a continuously-stirred tank reactor
%
% Syntax:  [problem,guess] = RayHicksCSTR
%
% Outputs:
%    problem - Structure with information on the optimal control problem
%    guess   - Guess for state, control and multipliers.
%
% Other m-files required: none
% Subfunctions: L (stageCost), 
%		E (boundaryCost), 
%		f (ODE right-hand side), 
%		g (path constraints), 
%		b (boundary constraints)
% MAT-files required: none
%
% Copyright (C) 2010 Paola Falugi, Eric Kerrigan and Eugene van Wyk. All Rights Reserved.
% This code is published under the BSD License.
% Department of Electrical and Electronic Engineering,
% Imperial College London London  England, UK 
% ICLOCS (Imperial College London Optimal Control) 5 May 2010
% iclocs@imperial.ac.uk

%------------- BEGIN CODE --------------

%Initial time. t0<tf
problem.time.t0=0;

% Final time. Let tf_min=tf_max if tf is fixed.
problem.time.tf_min=10;     
problem.time.tf_max=inf; 
guess.tf=120;

% Parameters bounds. pl=< p <=pu
problem.parameters.pl=[];
problem.parameters.pu=[];
guess.parameters=[];

% Initial conditions for system.
problem.states.x0=[0.9831 0.3918]; 

% Initial conditions for system. Bounds if x0 is free s.t. x0l=< x0 <=x0u
problem.states.x0l=problem.states.x0; 
problem.states.x0u=problem.states.x0l;

% State bounds. xl=< x <=xu
problem.states.xl=[0 0];
problem.states.xu=[1 1];

% Terminal state bounds. xfl=< xf <=xfu
problem.states.xfl=[0.2632 0.6519]; 
problem.states.xfu=[0.2632 0.6519];

% Guess the state trajectories with [x0 xf]
guess.states(:,1)=[0.9831 0.2632];
guess.states(:,2)=[0.3918 0.6519];

% Number of control actions N 
% Set problem.inputs.N=0 if N is equal to the number of integration steps.  
% Note that the number of integration steps defined in settings.m has to be divisible 
% by the  number of control actions N whenever it is not zero.

problem.inputs.N=40;       
      
% Input bounds
problem.inputs.ul=[0.0];
problem.inputs.uu=[2];

% Guess the input sequences with [u0 uf]
guess.inputs(:,1)=[0.0 455/600];


% Bounds for path constraint function gl =< g(x,u,p,t) =< gu
problem.constraints.gl=[];
problem.constraints.gu=[];

% Bounds for boundary constraints bl =< b(x0,xf,u0,uf,p,t0,tf) =< bu
problem.constraints.bl=[];
problem.constraints.bu=[];


% Choose the set-points if required
problem.setpoints.states=[0.2632 0.6519];
problem.setpoints.inputs=[455/600];

% store the necessary problem parameters used in the functions
problem.data=[];

% Get function handles and return to Main.m
problem.functions={@L,@E,@f,@g,@b};

%------------- END OF CODE --------------

function stageCost=L(x,xr,u,ur,p,t,data)

% L - Returns the stage cost.
% The function must be vectorized and
% xi, ui are column vectors taken as x(:,i) and u(:,i) (i denotes the i-th
% variable)
% 
% Syntax:  stageCost = L(x,xr,u,ur,p,t,data)
%
% Inputs:
%    x  - state vector
%    xr - state reference
%    u  - input
%    ur - input reference
%    p  - parameter
%    t  - time
%    data- structured variable containing the values of additional data used inside
%          the function
%
% Output:
%    stageCost - Scalar or vectorized stage cost
%
%  Remark: If the stagecost does not depend on variables it is necessary to multiply
%          the assigned value by t in order to have right vector dimesion when called for the optimization. 
%          Example: stageCost = 0*t;

%------------- BEGIN CODE --------------


c=x(:,1);T=x(:,2);u=u(:,1);
cr=xr(:,1);Tr=xr(:,2);

stageCost = 0.5*((c-cr).*(c-cr)+(T-Tr).*(T-Tr))+0.5*(u-ur).*(u-ur);

%------------- END OF CODE --------------


function boundaryCost=E(x0,xf,u0,uf,p,tf,data) 

% E - Returns the boundary value cost
%
% Syntax:  boundaryCost=E(x0,xf,u0,uf,p,tf,data)
%
% Inputs:
%    x0  - state at t=0
%    xf  - state at t=tf
%    u0  - input at t=0
%    uf  - input at t=tf
%    p   - parameter
%    tf  - final time
%    data- structured variable containing the values of additional data used inside
%          the function
%
% Output:
%    boundaryCost - Scalar boundary cost
%
%------------- BEGIN CODE --------------

boundaryCost=0;

%------------- END OF CODE --------------


function dx = f(x,u,p,t,data)

% f - Returns the ODE right hand side where x'= f(x,u,p,t)
% The function must be vectorized and
% xi, ui, pi are column vectors taken as x(:,i), u(:,i) and p(:,i). Each
% state corresponds to one column of dx.
% 
% 
% Syntax:  dx = f(x,u,p,t,data)
%
% Inputs:
%    x  - state vector
%    u  - input
%    p  - parameter
%    t  - time
%    data-structured variable containing the values of additional data used inside
%          the function 
%
% Output:
%    dx - time derivative of x
%
%  Remark: If the i-th ODE right hand side does not depend on variables it is necessary to multiply
%          the assigned value by a vector of ones with the same length  of t  in order 
%          to have  a vector with the right dimesion  when called for the optimization. 
%          Example: dx(:,i)= 0*ones(size(t,1)); 
%
%------------- BEGIN CODE --------------

%Ray-Hicks Model
% a=0.000195*600;q=10;En=25.2;k=300;Tc=2.9;Tf=3;

% Biegler's Model
a=0.000195*600;q=20;En=5;k=300;Tc=0.38158;Tf=0.3947;

c=x(:,1);T=x(:,2);u=u(:,1);

dx(:,1)=(1-c)/q-k*c.*exp(-En./T);
dx(:,2)=(Tf-T)/q+k*c.*exp(-En./T)-a*u.*(T-Tc);


%------------- END OF CODE --------------


function c=g(x,u,p,t,data)

% g - Returns the path constraint function where gl =< g(x,u,p,t) =< gu
% The function must be vectorized and
% xi, ui, pi are column vectors taken as x(:,i), u(:,i) and p(:,i). Each
% constraint corresponds to one column of c
% 
% Syntax:  c=g(x,u,p,t,data)
%
% Inputs:
%    x  - state vector
%    u  - input
%    p  - parameter
%    t  - time
%   data- structured variable containing the values of additional data used inside
%          the function
%
% Output:
%    c - constraint function
%
%------------- BEGIN CODE --------------

c=[];

%------------- END OF CODE --------------

function bc=b(x0,xf,u0,uf,p,tf,data)

% b - Returns a column vector containing the evaluation of the boundary constraints: bl =< bf(x0,xf,u0,uf,p,t0,tf) =< bu
%
% Syntax:  bc=b(x0,xf,u0,uf,p,tf,data)
%
% Inputs:
%    x0  - state at t=0
%    xf  - state at t=tf
%    u0  - input at t=0
%    uf  - input at t=tf
%    p   - parameter
%    tf  - final time
%    data- structured variable containing the values of additional data used inside
%          the function
%
%          
% Output:
%    bc - column vector containing the evaluation of the boundary function 
%
%------------- BEGIN CODE --------------

bc=[];

%------------- END OF CODE --------------






