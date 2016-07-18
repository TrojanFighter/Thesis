function res = bc(ya,yb,T)
alpha=0.65;
res = [ ya(1) - 4;
    ya(2) - 10;
    ya(3)-0.9;
    yb(2) - 2;
    yb(4);
    yb(6);
    alpha^2 + 2*(alpha +cos(yb(3))) * yb(5) - 1];
