%MAINMPC - Main script to solve the Model Predictive Control Problem
%
% Copyright (C) 2010 Paola Falugi, Eric Kerrigan and Eugene van Wyk. All Rights Reserved.
% This code is published under the BSD License.
% Department of Electrical and Electronic Engineering,
% Imperial College London London  England, UK 
% ICLOCS (Imperial College London Optimal Control) 5 May 2010
% iclocs@imperial.ac.uk


%--------------------------------------------------------------------

clear all 
format compact


[problem,guess]= testProblem;           % Fetch the problem definition
options= settings;                      % Get options and solver settings         
plant=@testPlant;                       % Get function handle of plant model

[infoNLP,data]=transcribeOCP(problem,guess,options); % Format for NLP solver
[nt,np,n,m,ng,nb,M,N,ns]=deal(data.sizes{:});
time=[];states=[];inputs=[];

P=20;                                   % Number of MPC iterations
      
            
% --- Begin MPC loop ---
for i=1:P
disp('Compute Control Action');disp(i);
    
solution = solveNLP(infoNLP,data);               % Solve the NLP


tc=solution.tf/N;                                % Control horizon
disp('Apply Control')
[x0,tv,xv,uv]=applyControl(tc,solution,plant,data,i); % Apply control 

time=[time;tv];                                 % Store results
states=[states;xv]; 
inputs=[inputs;uv];              

infoNLP.zl(nt+np+1:nt+np+n)=x0;                 % Update initial condition  
infoNLP.zu(nt+np+1:nt+np+n)=x0;  

infoNLP.z0=solution.z;                          % Update initial guess
data.x0t=x0.';
end
% --- End MPC loop ---


% Plot the solutions
figure(1)
hold on
plot(time,states)
title('States vs. time');
xlabel('Time')


figure(2)
hold on
plot(time,inputs)
title('Inputs vs. Time');
xlabel('Time')


figure(3)
hold on
plot(states(:,1),states(:,2))
title('Optimal state trajectory');
xlabel('x_1')
ylabel('x_2')

