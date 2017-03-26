function mi = plot(mi)
global dt

if ~mi.hit
for i = 1:size(mi.trail,1)
   if  isnan(mi.trail(i,3))
       mi.trail(i,3)= line(mi.position(1),mi.position(2),1,'Color',[1 0 0],'Marker','.','MarkerSize',1);
   else
       C = get(mi.trail(i,3),'Color');
       if C(2)+dt< 1
           C(2) = C(2)+dt;
           C(3) = C(3)+dt;
       elseif C(2)+dt/10 >= 1
           C(2) = 1;
           C(3) = 1;
       end
       set(mi.trail(i,3),'Color',C,'Marker','v','MarkerSize',2)%this was getting messed up by explosions;
   end
end

if mi.bodyaxis(1,2) >= 0
    headingangle = acos(mi.bodyaxis(1,1)/sqrt(mi.bodyaxis(1,1)^2+mi.bodyaxis(1,2)^2));
else
    headingangle = 2*pi-acos(mi.bodyaxis(1,1)/sqrt(mi.bodyaxis(1,1)^2+mi.bodyaxis(1,2)^2));
end

rotmatrix = [cos(headingangle) sin(headingangle); -sin(headingangle) cos(headingangle)];
delete(mi.graphic)
fv = mi.patchobj;
fv.vertices = [rotmatrix'*fv.vertices']';
fv.vertices(:,1) = fv.vertices(:,1) + mi.position(1);
fv.vertices(:,2) = fv.vertices(:,2) + mi.position(2);
mi.graphic = patch(fv);
set(mi.graphic,'EdgeColor',[0 0 0])
set(mi.graphic,'FaceColor','flat');


    
end

end