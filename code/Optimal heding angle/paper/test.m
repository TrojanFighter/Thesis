solinit = bvpinit(linspace(0,1),[2;3;1;2;1;2],6);

sol = bvp4c(@ode, @bc, solinit);
y = sol.y;
time = y(7)*sol.x;
%ut = -y(4,:);
plot(time,y(1,:),time,y(2,:),time,y(3,:)); %states vs time