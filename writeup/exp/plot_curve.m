% ====================================================
clc; clear all; close all; format short g; tic;
fprintf('Program Starts\n')
fprintf('==============\n')
% ====================================================
g = 1;
xA = 4;

alpha = [0.1:0.1:0.9];


% 1. Circle Data
% ==============
a = ((1 + g^2)/(1 - g^2)) * xA;
b = 0;
rA = ((2 * g)/abs(1 - g^2)) * xA;

N = 500; lb = -50; ub = 50;
theta = 0:pi/N:2*pi;
x_circle = rA * cos(theta) + a;
y_circle = rA * sin(theta) + b;

% 2. Polynomial Data
% ==================
xt = linspace(lb,ub,N);
yt = linspace(lb,ub,N);
[XT,YT] = meshgrid(xt,yt);
%legend('polynomial','AD Apollonius circle')
h = figure('position',[10 10 500 500]);


for i=1:length(alpha)
    
    Z = 16 * (g^2) * (alpha(i)^2)  * (xA^4) ...
        + 16 * (g^2) * (alpha(i)^2)  * (xA^2) * (XT.^2) ...
        - 32 * (g^2) * (alpha(i)^2)  * (xA^3) * (XT) ...
        + 16 * (g^2) * (alpha(i)^2)  * (xA^2) * (YT.^2) ...
        - (1-g^2)^2 * (xA^4 + (XT.^4) + (YT.^4) + 2*(xA^2)*(XT.^2) + 2*(xA^2)*(YT.^2) + 2*(XT.^2).*(YT.^2)) ...
        - (1+g^2)^2 * 4 *(xA^2) * (XT.^2) ...
        - 16 * (g^4) * (alpha(i)^4)  * (xA^4) ...
        + 2 * (1 - g^4) * (2*(xA^3)*(XT) + 2*(xA)*(XT.^3) + 2*(xA)*(XT).*(YT.^2)) ...
        - 8 * (g^2) * (1-g^2) * (alpha(i)^2) * (xA^4 + ((xA^2) * (XT.^2)) + ((xA^2) * (YT.^2))) ...
        +16 * (g^2) * (1+g^2) * (alpha(i)^2) * (xA^3) * XT;
    % =========================================
    % Plot
    
    % polynomial
    v = [0,0];
    
    %[C,h2]= contour(XT,YT,Z,v,'linewidth',3,'ShowText','on');
    C = contourcs(xt,yt,Z,v);
    %plot(C(1).X,C(1).Y )
    %contour3(XT,YT,Z);
    % hold on
    % meshc(XT,YT,Z)
    plot(C(2).X,C(2).Y)
    hold on
    
end


% Circle
plot(x_circle,y_circle,'linewidth',3); %hold on
colormap jet
xlabel('x_T'); ylabel('y_T'); zlabel('p(x_T,y_T)')
set(gcf,'color','w'); set(gca,'FontSize',10);
box on; axis square; grid on
axis ([lb ub lb ub]);

legend('\alpha = 0','\alpha = 0.1','','','','','','','\alpha = 1','Location','northwest')
%h = figure('position',[10 10 500 500]);

% =========================================
saveTightFigure(h,'polynomial+circle.pdf')
% =========================================
fprintf('==============\n')
fprintf('Program Ends\n')
toc
% =========================================