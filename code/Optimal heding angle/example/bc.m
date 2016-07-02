function res = bc(ya,yb,T)
res = [ ya(1) - 1; ya(2) - 2; yb(1) - 3; yb(4);
yb(3)*yb(2)-0.5*yb(4)^2];