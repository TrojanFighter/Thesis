% ODE’s of augmented states
function dydt = ode(t,y,T)
dydt = T*[ y(2);-y(4);0;-y(3)];
