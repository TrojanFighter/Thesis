function dx = rocket(t,x)

alpha = 0.65;
%% Defining the optimal heading angles to be used in the simulation
% psi = asin(x(6)/(x(2)*sqrt(x(5)^2 + (x(6)/x(2))^2)));
% phi = asin(x(6)/(x(1)*sqrt((1-x(4))^2 + (x(6)/x(1))^2)));
psi = atan2(x(6),x(2) * x(5));
phi = atan2(x(6),x(1)*(1-x(4)));
Xs = (1-x(4))*sin(x(3)) - x(6)/x(1) * cos(x(3)) + x(6)/x(2);
Xc = (1-x(4))*cos(x(3)) + x(6)/x(1) * sin(x(3)) - x(5);
% xi = asin(Xs/sqrt(Xs^2 + Xc^2));
xi = atan2(Xs,Xc);
%% The change in system dynamics (the x_dot vector)
dx(1) = alpha*cos(phi) - cos(x(3)-xi);
dx(2) = -cos(xi)-cos(psi);
dx(3) = -alpha/x(1) * sin(phi) + 1/x(1) * sin(x(3)-xi)- 1/x(2)*sin(psi) + 1/x(2) * sin(xi);
dx(4) = x(6)/((x(1)^2))*(sin(x(3)-xi)-alpha*sin(phi));
dx(5) = x(6)/((x(2)^2))*(sin(xi)-sin(psi));
dx(6) = (1-x(4))*sin(x(3)-xi) - x(6)/x(1)*cos(x(3)-xi);

dx = dx';