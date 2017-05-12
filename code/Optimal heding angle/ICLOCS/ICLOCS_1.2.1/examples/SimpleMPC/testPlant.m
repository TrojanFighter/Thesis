function dx = testPlant(x,u,p,t)
% testPlant - Returns ODE rhs of the PLANT where x'= f(x,u,p,t)
% For direct transcription methods: the function must be vectorized and
% xi, ui, pi are column vectors taken as x(:,i), u(:,i) and p(:,i). Each
% state corresponds to one column of dx.
% 
% 
% Syntax:  dx = testPlant(x,u,p,t)
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
% Copyright (C) 2010 Paola Falugi, Eric Kerrigan and Eugene van Wyk. All Rights Reserved.
% This code is published under the BSD License.
% Department of Electrical and Electronic Engineering,
% Imperial College London London  England, UK 
% ICLOCS (Imperial College London Optimal Control) 5 May 2010
% iclocs@imperial.ac.uk

%------------- BEGIN CODE --------------

x1=x(:,1);
x2=x(:,2);

u=u(:,1);

dx(:,1)=1.2*x2+0.1*sin(t);
dx(:,2)=0.2*sin(x1)+u;
%------------- END OF CODE --------------