% MAIN - Main script to solve a Standard MPC Problem for discrete time systems
%
% Copyright (C) 2010 Paola Falugi, Eric Kerrigan and Eugene van Wyk. All Rights Reserved.
% This code is published under the BSD License.
% Department of Electrical and Electronic Engineering,
% Imperial College London London  England, UK 
% ICLOCS (Imperial College London Optimal Control) 5 May 2010
% iclocs@imperial.ac.uk

%--------------------------------------------------------------------

clear functions;clear all


format compact


%color options for the pictures
colMPC='b--';

P1=100;                                   % Number of iterations 


[pb,guess]=Discrete_Sys;       % Fetch the problem definition
opts= settings_Dis;          % Get options and solver settings


[infoNLP,data_mpc]=transcribeOCP(pb,guess,opts);% Format for NLP solver

% Get some dimensions from the problem definitions

[nt,np,n,m,ng,nb,M,N]=deal(data_mpc.sizes{1:8});



load termset

% -----------------------------
 
time=[];states=[]; inputs=[];cost=[];status_mpc=[];

                                     % Number of iterations  
    
    for i=1:P1
    disp('Standard MPC computation. Iteration:');disp(i);

    [sol_mpc,status] = solveNLP(infoNLP,data_mpc);  % Solve the NLP for reference

    
    
    y0=sol_mpc.X((M-1)/N+1,:);                   % Next y0

    infoNLP.zl(nt+np+1:nt+np+n)=y0;        % Update initial condition  for reference
    infoNLP.zu(nt+np+1:nt+np+n)=y0;  

    infoNLP.z0=sol_mpc.z;                  % Update initial guess for reference
    
    
    states=[states;sol_mpc.X(1:(M-1)/N,:)];
    inputs=[inputs;sol_mpc.U(1:(M-1)/N,:)]; 
    time=[time;(i-1)];
    cost=[cost;costFunction(sol_mpc.z,data_mpc)];
    status_mpc=[status_mpc;status];
    data_mpc.x0t=y0.';

    
    end
 
% Plot Solutions

figure(1)
 subplot(1,2,1)
 plot(time,states(:,1),colMPC);
 hold on
 title('MPC');ylabel('x_1');
 xlabel('Time')
 
 
subplot(1,2,2)
 plot(time,states(:,2),colMPC);
 hold on
 title('MPC');ylabel('x_2');
 xlabel('Time')

figure(2)
plot(time,inputs(:,1),colMPC);
hold on
title('MPC');ylabel('Input');
xlabel('Time')

figure(3)
plot(cost);
hold on
title('Cost')
xlabel('Steps')

figure(4)
plot(states(:,1),states(:,2),colMPC);  
hold on
xlabel('x_1')
ylabel('x_2')

 
 
    

