function dydt = ode(t,y,T)

alpha=0.65;
psi = atan2(y(6),y(2) * y(5));
phi = atan2(y(6),y(1)*(1-y(4)));
Xs = (1-y(4))*sin(y(3)) - y(6)/y(1) * cos(y(3)) + y(6)/y(2);
Xc = (1-y(4))*cos(y(3)) + y(6)/y(1) * sin(y(3)) - y(5);
xi = atan2(Xs,Xc);


dydt = T*[ alpha*cos(phi) - cos(y(3))-xi;
-cos(xi)-cos(psi);
 -alpha/y(1) * sin(phi) + 1/y(1) * sin(y(3)-xi)- 1/y(2)*sin(psi) + 1/y(2) * sin(xi);
 y(6)/((y(1)^2))*(sin(y(3)-xi)-alpha*sin(phi));
 y(6)/((y(2)^2))*(sin(xi)-sin(psi));
 (1-y(4))*sin(y(3)-xi) - y(6)/y(1)*cos(y(3)-xi)];