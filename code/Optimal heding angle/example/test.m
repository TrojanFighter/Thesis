solinit = bvpinit(linspace(0,1),[2;3;1;1],2);
sol = bvp4c(@ode, @bc, solinit);
y = sol.y;
time = y(5)*sol.x;
ut = -y(4,:);
plot(time,y(1,:),time,y(2,:));