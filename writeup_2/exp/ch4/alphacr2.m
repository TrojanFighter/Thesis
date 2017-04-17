% ====================================================
clc; clear all; close all; format short g; tic;
fprintf('Program Starts\n')
fprintf('==============\n')
% ====================================================
gamma=1.25;
a=[-3:1:3]; % x_T/x_A
b=[-3:1:3];  % y_T/x_A

for i = 1: length(b)
    % calculate critical alpha:
    alpha_cr(i,:) =(1/(2*gamma))*(- sqrt((1-a).^2+b(i)^2)+ (gamma* sqrt((1+a).^2+b(i)^2)));
    
    % find elements less than zero or larger than one.
    k = find ((alpha_cr(i,:)<0)|(alpha_cr(i,:)>1));
    % set these elements to zero.
    alpha_cr (i,k)= 0;
    
    
    
    %=================
    
    xlabel('x_T/x_A');
    ylabel('\alpha_c_r');
    title('critical normalized speed for \gamma=1')
    plot (a, alpha_cr(i,:))
    grid on
    hold on
end


% =========================================
fprintf('==============\n')
fprintf('Program Ends\n')
toc
% =========================================