clc; close all; clear; format short g

solinit = bvpinit(linspace(0,1,5),[4;10;0.9;0.45;0.25;0.25],6);

%options = bvpset('RelTol',0.1);
%options = bvpset('NMax', 1E4);
sol = bvp4c(@ode, @bc, solinit);

y = sol.y;
x = y(7) * sol.x;

%ut = -y(4,:)

% ===========
h = figure('position',[100, 100, 500, 750]); set(gcf,'color','w'); set(gca,'FontSize',24);
plot(x,y(1,:),'-b',x,y(2,:),'-k',x,y(3,:),'-r','LineWidth',3); 
grid on
h_legend = legend('[y_1]: R','[y_2]: r','[y_3]: \theta','Location','NorthEast');
set(h_legend,'FontSize',14); set(gca,'FontSize',14);

xlabel('Time [sec]','FontSize', 16); ylabel('States','FontSize', 16)
saveTightFigure(h,'temp.png')