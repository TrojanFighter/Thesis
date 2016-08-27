clc; close all; clear; format short g
% ===========
%
% Code to replicate figures from:
% "Active Target Defense Differential Game" - Pachter, Garcia, Casbeer 
% Allerton 2014
%
% ===========
solinit = bvpinit(linspace(0,1.75,3),[4;10;0.9;0.45;0.25;0.25],6);
sol = bvp4c(@ode, @bc, solinit);

% ===========
% Figure 7a 

y = sol.y;
x = y(7) * sol.x;

h = figure('position',[100, 100, 500, 750]); set(gcf,'color','w'); set(gca,'FontSize',24);
plot(x,y(1,:),'-b',x,y(2,:),'-k',x,y(3,:),'-r','LineWidth',3); 
grid on
h_legend = legend('[y_1]: R','[y_2]: r','[y_3]: \theta','Location','NorthEast');
set(h_legend,'FontSize',14); set(gca,'FontSize',14);

xlabel('Time [sec]','FontSize', 16); ylabel('States','FontSize', 16)
saveTightFigure(h,'temp.png')

% ===========
% Figure 7b 

psi = atan2(y(6,:),y(2,:) .* y(5,:));
phi = atan2(y(6,:),y(1,:) .* (1-y(4,:)));
Xs = (1-y(4,:)) .* sin(y(3,:)) - y(6,:)./y(1,:) .* cos(y(3,:)) + y(6,:)./y(2,:);
Xc = (1-y(4,:)) .* cos(y(3,:)) + y(6,:)./y(1,:) .* sin(y(3,:)) - y(5,:);
xi = atan2(Xs,Xc);

h = figure('position',[100, 100, 500, 750]); set(gcf,'color','w'); set(gca,'FontSize',24);
plot(x,phi,'-b',x,psi,'-k',x,xi,'-r','LineWidth',3); 
grid on
h_legend = legend('\phi^*','\psi^*','\chi^*','Location','NorthWest');
set(h_legend,'FontSize',14); set(gca,'FontSize',14);

xlabel('Time [sec]','FontSize', 16); ylabel('States','FontSize', 16)
saveTightFigure(h,'temp.png')

% =========== % Figure 8 
alpha = 0.65; beta = 1;
VA = 1; VT = alpha * VA; VD = beta * VA;

xA(1) = 5;     yA(1) = 0;
xD(1) = -5;    yD(1) = 0;
xT(1) = 2.514; yT(1) = 3.133;
lambda(1) = pi + (atan((yT(1)-yA(1))/(xT(1)- xA(1))));

phi_hat(1) = phi(1) + lambda(1);
psi_hat(1) = psi(1) + y(3,1) + lambda(1) - pi;
xi_hat(1) =  lambda(1) + y(3,1) - xi(1);
%phi_hat(1) = 2.7033; psi_hat(1) = 0.7115; xi_hat(1)  = 2.4301;

xTdot(1) = VT * cos(phi_hat(1)); yTdot(1) = VT * sin(phi_hat(1));
xAdot(1) = VA * cos(xi_hat(1));  yAdot(1) = VA * sin(xi_hat(1));
xDdot(1) = VD * cos(psi_hat(1)); yDdot(1) = VD * sin(psi_hat(1));

for i = 1 : length(x) - 1,
    xT(i + 1) = xT(i) + xTdot(i); yT(i + 1) = yT(i) + yTdot(i);
    xA(i + 1) = xA(i) + xAdot(i); yA(i + 1) = yA(i) + yAdot(i);
    xD(i + 1) = xD(i) + xDdot(i); yD(i + 1) = yD(i) + yDdot(i);
    
    lambda(i + 1) = pi + atan((yT(i+1)-yA(i+1))/(xT(i+1)- xA(i+1)));
    
    phi_hat(i+1) = phi(i+1) + lambda(i+1);
    psi_hat(i+1) = psi(i+1) + y(3,i+1) + lambda(i+1) - pi;
    xi_hat(i+1) =  lambda(i) + y(3,i+1) - xi(i+1);
    %phi_hat(i+1) = 2.7033; psi_hat(i+1) = 0.7115; xi_hat(i+1)  = 2.4301;

    xTdot(i+1) = VT * cos(phi_hat(i+1)); yTdot(i+1) = VT * sin(phi_hat(i+1));
    xAdot(i+1) = VA * cos(xi_hat(i+1));  yAdot(i+1) = VA * sin(xi_hat(i+1));
    xDdot(i+1) = VD * cos(psi_hat(i+1)); yDdot(i+1) = VD * sin(psi_hat(i+1));
end

h = figure('position',[100, 100, 1000, 750]); set(gcf,'color','w'); set(gca,'FontSize',24);
plot(xT(1),yT(1),'-ob','MarkerSize',12,'MarkerFaceColor','b'); hold on
plot(xA(1),yA(1),'-or','MarkerSize',12,'MarkerFaceColor','r'); 
plot(xD(1),yD(1),'-ok','MarkerSize',12,'MarkerFaceColor','k'); 
plot(xT,yT,'-b',xA,yA,'-r',xD,yD,'-k','LineWidth',3); 
grid on
h_legend = legend('Target','Attacker','Defender','Location','North');
set(h_legend,'FontSize',14); set(gca,'FontSize',14);
xlabel('x','FontSize', 16); ylabel('y','FontSize', 16)
%axis([-8 8 -2 8])
saveTightFigure(h,'temp.png')
% ===========
