function [solution,status]=solveNLP(NLP,data)

%SOLVENLP - Solve the nonlinear program and return the solution
%
% Syntax:  solution = solveNLP(NLP,options)
%
% Inputs:
%    NLP     - Information required by the NLP solver
%    data    - Data passed to the functions evaluated during optimization
%
% Outputs:
%    solution - Structure containing the optimal final time, parameters,
%               states, controls and Lagrange multipliers.
%    status - variable containing the exit condition of the solver. See the
%             fmincon and ipopt description 
%
% Other m-files required: ipopt.mex and or fmincon.m
% Subfunctions: none
% MAT-files required: none
%
% Copyright (C) 2010 Paola Falugi, Eric Kerrigan and Eugene van Wyk. All Rights Reserved.
% This code is published under the BSD License.
% Department of Electrical and Electronic Engineering,
% Imperial College London London  England, UK 
% ICLOCS (Imperial College London Optimal Control) 5 May 2010
% iclocs@imperial.ac.uk

%------------- BEGIN CODE --------------

% Solve the nonlinear program:


[nt,np,n,m,ng,nb,M,N,ns]=deal(data.sizes{:});
NLP.z0(nt+np+1:n+nt+np)=data.x0t;


switch(data.options.NLPsolver)

    case{'ipopt'} 
       opt=data.multipliers;              % Solve the NLP using IPOPT
       opt.ipopt=data.options.ipopt;
       opt.cl=NLP.cl;                    % Lower bounds on constraints.
       opt.cu=NLP.cu;                      % Upper bounds on constraints
       opt.lb=NLP.zl;                       % Lower bound on the variables.
       opt.ub=NLP.zu;                      % Upper bound on the variables.
       opt.auxdata=data;
       
       tA=cputime;
       [z, info] = ipopt_auxdata(NLP.z0,data.funcs,opt);
       tB=cputime;
       solution.multipliers.zl=info.zl;
       solution.multipliers.zu=info.zu;
       solution.multipliers.lambda=info.lambda;
       status=info.status;
    
    case{'fmincon'}                               % Solve the NLP using FMINCON
        
        tA=cputime;
        [z,cost,status,output,multipliers]=fmincon(@(z)fminCost(z,data),NLP.z0,...
            [],[],[],[],NLP.zl,NLP.zu,@(z)fminConst(z,data,NLP),data.options.fmincon);
        tB=cputime;
        ni=output.iterations;
        solution.multipliers=multipliers;
        solution.iterates=ni;
    otherwise
        disp('Unknown NLP solver. Check spelling.');

end


% Store the results in solution structure:
%solution.status=status;

solution.computation_time=tB-tA;
solution.z=z;

if nt;solution.tf=z(1); else solution.tf=data.tf; end
if np;solution.p=z(nt+1:nt+np); else solution.p=[]; end
t0=data.t0;

if strcmp(data.options.transcription,'multiple_shooting')


    solution.X=reshape(data.map.Vx*z,n,M)';
    solution.x0=solution.X(1,:)';
    usp=data.map.Vu*z;
    solution.U=reshape([usp;usp(end-m+1:end)],m,M)';
    solution.T=(solution.tf-t0)*[0;cumsum(data.tau)]+t0; 

else
    
    solution.X=reshape(data.map.Vx*z,n,M)';
    solution.x0=solution.X(1,:)';
    usp=reshape(data.map.Vu*z,m,N)';
    solution.U=kron(usp,ones((M-1)/N,1));
    solution.U=[solution.U;solution.U(end,:)];
    solution.T=(solution.tf-t0)*[0;cumsum(data.tau)*data.Nm/ns]+data.k0;
     
  
    if strcmp(data.options.transcription,'hermite')
        solution.X=solution.X(1:2:end,:);
        solution.U=solution.U(1:2:end,:);
        solution.T=solution.T(1:2:end);

    end


end

%------------- END OF CODE --------------