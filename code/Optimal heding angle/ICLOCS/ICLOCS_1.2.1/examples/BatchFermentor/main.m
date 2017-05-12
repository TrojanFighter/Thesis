% MAIN - Main script to solve the Optimal Control Problem
%
% Copyright (C) 2010 Paola Falugi, Eric Kerrigan and Eugene van Wyk. All Rights Reserved.
% This code is published under the BSD License.
% Department of Electrical and Electronic Engineering,
% Imperial College London London  England, UK 
% ICLOCS (Imperial College London Optimal Control) 5 May 2010
% iclocs@imperial.ac.uk

%--------------------------------------------------------

clear all;format compact;

global sol;               % Initialize solution structure used in IPOPT
[problem,guess]=BatchFermentor;       % Fetch the problem definition
options= settings;                      % Get options and solver settings         

[infoNLP,data]=transcribeOCP(problem,guess,options); % Format for NLP solver

% Initialize the dual point
nc=size(data.jacStruct,1);
[nt,np,n,m,M,N]=deal(data.sizes{[1:4,7:8]});
nz=nt+np+n*M+m*N;
data.multipliers.zl=2*ones(1,nz);
data.multipliers.zu=2*ones(1,nz);
data.multipliers.lambda=2*ones(1,nc);
 

[solution,status] = solveNLP(infoNLP,data);      % Solve the NLP

output(solution,options,data);         % Output solutions

