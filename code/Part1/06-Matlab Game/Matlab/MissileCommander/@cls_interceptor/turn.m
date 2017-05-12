function mi = turn(mi)
bodyaxis = mi.bodyaxis;
bodyhandle = mi.graphic;
position = mi.position;
[headingangle, angle, yb, rotmatrix] = getheading(mi);
%turn the body


% rotmatrix = [cos(theta) sin(theta); -sin(theta) cos(theta)];
% bodyaxis = [rotmatrix*bodyaxis]';
%     if yb >= 0 
%         theta = 3*pi/180;
%     else
%         theta = -3*pi/180;
%     end
if angle < 20*pi/180
   mi.thrusting = true;
else
   mi.thrusting = false;
end

if yb >= 0
    dir = 1;
else
    dir = -1;
end

if angle > 10*pi/180
    theta = dir*10*pi/180;
else
    theta = dir*angle;
end

headingangle = headingangle+theta; 
bodyaxis(1,1) = cos(headingangle);
bodyaxis(1,2) = sin(headingangle);
bodyaxis(2,1) = sin(headingangle);
bodyaxis(2,2) = cos(headingangle);



% delete(bodyhandle)
% fv = mi.patchobj;
% fv.vertices = [rotmatrix'*fv.vertices']';
% fv.vertices(:,1) = fv.vertices(:,1) + position(1);
% fv.vertices(:,2) = fv.vertices(:,2) + position(2);
% 
% mi.graphic = patch(fv);
% set(mi.graphic,'EdgeColor',[.5 .5 .5])
% set(mi.graphic,'FaceColor','flat');


mi.bodyaxis = bodyaxis;







end