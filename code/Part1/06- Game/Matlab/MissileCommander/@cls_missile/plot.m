function m = plot(m)
global dt

if ~m.hit
for i = 1:size(m.trail,1)
   if  isnan(m.trail(i,3))
       m.trail(i,3)= line(m.position(1),m.position(2),1,'Color',[1 0 0],'Marker','.','MarkerSize',2);
   else
       C = get(m.trail(i,3),'Color');
       if C(2)+dt< 1
           C(2) = C(2)+dt;
           C(3) = C(3)+dt;
       elseif C(2)+dt/10 >= 1
           C(2) = 1;
           C(3) = 1;
       end
       set(m.trail(i,3),'Color',C,'Marker','v','MarkerSize',2)%this was getting messed up by explosions;
   end
end

if m.velocity(1,2) >= 0
    headingangle = acos(m.velocity(1,1)/sqrt(m.velocity(1,1)^2+m.velocity(1,2)^2));
else
    headingangle = 2*pi-acos(m.velocity(1,1)/sqrt(m.velocity(1,1)^2+m.velocity(1,2)^2));
end

rotmatrix = [cos(headingangle) sin(headingangle); -sin(headingangle) cos(headingangle)];
delete(m.graphic)
fv = m.patchobj;
fv.vertices = [rotmatrix'*fv.vertices']';
fv.vertices(:,1) = fv.vertices(:,1) + m.position(1);
fv.vertices(:,2) = fv.vertices(:,2) + m.position(2);
m.graphic = patch(fv);
set(m.graphic,'EdgeColor',[0 0 0])
set(m.graphic,'FaceColor','flat');
    
end

end