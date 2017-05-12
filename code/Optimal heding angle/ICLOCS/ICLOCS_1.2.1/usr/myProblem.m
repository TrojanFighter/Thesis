function [problem,guess] = myProblem

%myProblem - Template file for optimal control problem definition
%
%Syntax:  [problem,guess] = myProblem
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


% Initial time. t0<tf.  For discrete time systems is the initial index
problem.time.t0=t0;


% Final time. Let tf_min=tf_max if tf is fixed.
problem.time.tf_min=final_time_min;     
problem.time.tf_max=final_time_max; 
guess.tf=final_time_guess;

% Parameters bounds. pl=< p <=pu
problem.parameters.pl=[p1_lowerbound ...];
problem.parameters.pu=[p1_upperbound ...];
guess.parameters=[p1_guess p2_guess ...];

% Initial conditions for system. Bounds if x0 is free s.t. x0l=< x0 <=x0u
problem.states.x0l=[x1(t0)_lowerbound ... xn(t0)_lowerbound]; 
problem.states.x0u=[x1(t0)_upperbound ... xn(t0)_upperbound]; 

% State bounds. xl=< x <=xu
problem.states.xl=[x1_lowerbound ... xn_lowerbound];
problem.states.xu=[x1_upperbound ... xn_upperbound];

% Terminal state bounds. xfl=< xf <=xfu
problem.states.xfl=[x1(tf)_lowerbound ... xn(tf)_lowerbound]; 
problem.states.xfu=[x1(tf)_upperbound ... xn(tf)_upperbound];

% Guess the state trajectories with [x0 xf]
guess.states(:,1)=[x1(t0) x1(tf)];
% ...
guess.states(:,n)=[xn(t0) xn(tf)];


% Number of control actions N 
% Set problem.inputs.N=0 if N is equal to the number of integration steps.  
% Note that the number of integration steps defined in settings.m has to be divisible 
% by the  number of control actions N whenever it is not zero.

problem.inputs.N=0;

% Input bounds
problem.inputs.ul=[u1_lowerbound ... um_lowerbound];
problem.inputs.uu=[u1_upperbound ... um_upperbound];

% Guess the input sequences with [u0 uf]
guess.inputs(:,1)=[u1(t0) u1(tf)]];
%...
guess.inputs(:,m)=[um(t0) um(tf)]];


% Choose the set-points if required
problem.setpoints.states=[x1_setpoint ... xn_setpoint];
problem.setpoints.inputs=[u1_setpoint ... um_setpoint];

% Bounds for path constraint function gl =< g(x,u,p,t) =< gu
problem.constraints.gl=[g1_lowerbound g2_lowerbound ...];
problem.constraints.gu=[g1_upperbound g2_upperbound ...];

% Bounds for boundary constraints bl =< b(x0,xf,u0,uf,p,t0,tf) =< bu
problem.constraints.bl=[b1_lowerbound b2_lowerbound ...];
problem.constraints.bu=[b1_upperbound b2_upperbound ...];


% Store the necessary problem parameters used in the functions
problem.data=[];


% Get function handles and return to Main.m
problem.functions={@L,@E,@f,@g,@b};

%------------- END OF CODE --------------

function stageCost=L(x,xr,u,ur,p,t)

% L - Returns the stage cost.
% The function must be vectorised and
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
%    stageCost - Scalar or vectorised stage cost
%
%  Remark: If the stagecost does not depend on variables it is necessary to multiply
%          the assigned value by t in order to have right vector dimesion when called for the optimisation. 
%          Example: stageCost = 0*t;

%------------- BEGIN CODE --------------


%Define states and setpoints
x1 = x(:,1); xr1=xr(:,1);
%...
xn=x(:,n); xrn=xr(:,n);

%Define inputs
u1 = u(:,1);ur1=ur(:,1);
% ...
um = u(:,m);urm=ur(:,m);

stageCost = L(x1,xr1...xn,u1,ur1...um,p,t);

%------------- END OF CODE --------------


function boundaryCost=E(x0,xf,u0,uf,p,tf) 


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


boundaryCost=E(x0,xf,u0,uf,p,tf);

%------------- END OF CODE --------------


function dx = f(x,u,p,t)

% f - Returns the ODE right hand side where x'= f(x,u,p,t)
% The function must be vectorised and
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
%          to have  a vector with the right dimesion  when called for the optimisation. 
%          Example: dx(:,i)= 0*ones(size(t,1)); 
%
%------------- BEGIN CODE --------------


%Define states
x1 = x(:,1);
%...
xn=x(:,n),

%Define inputs
u1 = u(:,1);
% ...
um = u(:,m);


%Define ODE right-hand side

dx(:,1) = f1(x1,..xn,u1,..um,p,t);
%...
dx(:,n) = fn(x1,..xn,u1,..um,p,t);

%------------- END OF CODE --------------


function c=g(x,u,p,t)

% g - Returns the path constraint function where gl =< g(x,u,p,t) =< gu
% The function must be vectorised and
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


%Define states
x1 = x(:,1);
%...
xn=x(:,n),

%Define inputs
u1 = u(:,1);
% ...
um = u(:,m);

c(:,1)=g1(x1,...,u1,...p,t);
c(:,2)=g2(x1,...,u1,...p,t);



%------------- END OF CODE --------------

function bc=b(x0,xf,u0,uf,p,tf)

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


bc(1,:)=b1(x0,xf,u0,uf,p,tf);
bc(2,:)=b2(x0,xf,u0,uf,p,tf);
%...
%------------- END OF CODE --------------






