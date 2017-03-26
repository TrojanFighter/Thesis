function [headingangle, theta, yb, rotmatrix] = getheading(mi)

bodyaxis = mi.bodyaxis;
position = mi.position;
target   = mi.target;
%v1 = [target(1) - position(1), target(2) - position(2)];
v1 = target - position;
%magv1 = sqrt(v1(1)^2+v1(2)^2);
v2 = [bodyaxis(1,1),bodyaxis(1,2)];
theta = acos(dot(v1,v2)/norm(v1));

if bodyaxis(1,2) >= 0
    headingangle = acos(bodyaxis(1,1));%(equals one) /sqrt(bodyaxis(1,1)^2+bodyaxis(1,2)^2));
else
    headingangle = 2*pi-acos(bodyaxis(1,1));%(equals one)/sqrt(bodyaxis(1,1)^2+bodyaxis(1,2)^2));
end

rotmatrix = [cos(headingangle) sin(headingangle); -sin(headingangle) cos(headingangle)];
%bftarget = [rotmatrix*[target(1)-position(1),target(2)- position(2)]']';
bftarget = [rotmatrix*v1']';

yb = bftarget(2);



end