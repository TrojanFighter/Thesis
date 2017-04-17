% ====================================================
clc; clear all; close all; format short g; tic;
fprintf('Program Starts\n')
fprintf('==============\n')
% ====================================================
gamma=1;
a=[-3:1:3]; % x_T/x_A
b= 1;  % y_T/x_A
%i=1;

for (i=1:length(a))
    
    alpha_cr(1,i)=(1/(2*gamma))*(- sqrt((1-a(i))^2+b^2)+ (gamma* sqrt((1+a(i))^2+b^2)));
    
    if (alpha_cr(1,i)<1)& (alpha_cr(1,i)>0)
       alpha_cr(1,i)=alpha_cr(1,i);
    else
        alpha_cr(1,i)=0;
    end
    
   % plot (a(i), alpha_cr(i))
end
plot (a, alpha_cr)
grid on
xlabel('x_T/x_A');
ylabel('\alpha_c_r');
title('critical normalized speed for \gamma=0.8')

% =========================================
fprintf('==============\n')
fprintf('Program Ends\n')
toc
% =========================================