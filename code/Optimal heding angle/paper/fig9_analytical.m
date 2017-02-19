clc; close all; clear; format short g
% ===========
%
% Code to replicate figures from:
% "Active Target Defense Differential Game" - Pachter, Garcia, Casbeer 
% Allerton 2014
%
% ===========
% Figure 8 

alpha = 0.415; 
beta = 1;

VA = 1; 
VT = alpha * VA; 
VD = beta * VA;

xA = 5;     yA = 0;
xD = -5;    yD = 0;
xT = 2.514; yT = 3.133;

phi_hat = 2.7033;
psi_hat = 0.7115;
xi_hat  = 2.4301;

%syms y
c4 = 1 - alpha^2;
c3 = -2*(1 - alpha^2) * yT;
c2 = (1 - alpha^2)*yT^2 + xA^2 - alpha^2*xT^2;
c1 = -2 * xA^2 * yT;
c0 = xA^2 * yT^2;

y = roots([c4 c3 c2 c1 c0])
%solve( c4 * y^4 + c3 * y^3 + c2 * y^2 + c1 * y + c0 == 0, y)

tf = (abs(xD)/cos(psi_hat)) / (VD)

alpha_bar = (sqrt((xA+xT)^2+(yT^2)) - sqrt((xA-xT)^2+(yT^2)))/(2*xA)

Target_distance = VT * tf

xT = [xT xT + Target_distance*cos(phi_hat)]; 
xA = [xA 0];
xD = [xD 0];

yT = [yT yT + Target_distance*sin(phi_hat)];
yA = [yA real(y(3))];
yD = [yD real(y(3))];

h = figure('position',[100, 100, 1000, 750]); set(gcf,'color','w'); set(gca,'FontSize',24);

plot(xT,yT,'-b','lineWidth',3); hold on
plot(xA,yA,'-r','lineWidth',3)
plot(xD,yD,'-k','lineWidth',3)

plot(xT(1),yT(1),'-ob','MarkerSize',12,'MarkerFaceColor','b'); hold on
plot(xA(1),yA(1),'-or','MarkerSize',12,'MarkerFaceColor','r'); 
plot(xD(1),yD(1),'-ok','MarkerSize',12,'MarkerFaceColor','k');

txt1 = 'T'; txt2 = 'A'; txt3 = 'D';
tx1 = text(xT(1)+0.2,yT(1)+0.2,txt1); tx1.FontSize = 20;
tx2 = text(xA(1)-0.2,yA(1)-0.2,txt2); tx2.FontSize = 20;
tx3 = text(xD(1)+0.2,yD(1)-0.2,txt3); tx3.FontSize = 20;

grid on
h_legend = legend('Target','Attacker','Defender','Location','NorthEast');
set(h_legend,'FontSize',14); set(gca,'FontSize',14);
xlabel('x','FontSize', 16); ylabel('y','FontSize', 16)
axis([-5 5 -1 5])
saveTightFigure(h,'fig9.png')
% ===========