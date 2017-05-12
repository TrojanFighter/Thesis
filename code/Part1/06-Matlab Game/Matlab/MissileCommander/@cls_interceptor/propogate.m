function mi = propogate(mi)
global dt
global round
global launchplayer;
if mi.launched
   mi.age = mi.age + dt;
   if mi.age == .1
       stop(launchplayer);
       play(launchplayer);
   end
   
   if ~mi.hit && ~mi.stop
       
       
       if mi.age > .1
       mi = turn(mi);
           if mi.thrusting
                thrust = mi.bodyaxis(1,:)*.085*((round+2)/3);
           else
                thrust = [0 0];
           end
       mi.velocity(1) = mi.velocity(1)+ thrust(1);
       mi.velocity(2) = mi.velocity(2)+ thrust(2)- .03*dt;
       mi.position(1) = mi.position(1)+mi.velocity(1)*dt;
       mi.position(2) = mi.position(2)+mi.velocity(2)*dt;
       else
         mi.velocity(2) = mi.velocity(2)+ .085 - .03*dt; 
         mi.position(2) = mi.position(2)+mi.velocity(2)*dt;
%         delete(mi.graphic)
%         fv = mi.patchobj;
%         rotmatrix = [cos(pi/2) sin(pi/2); -sin(pi/2) cos(pi/2)];
%         fv.vertices = [rotmatrix'*fv.vertices']';
%         fv.vertices(:,1) = fv.vertices(:,1) + mi.position(1);
%         fv.vertices(:,2) = fv.vertices(:,2) + mi.position(2);
%         mi.graphic = patch(fv);
%         set(mi.graphic,'EdgeColor',[.5 .5 .5])
%         set(mi.graphic,'FaceColor','flat');
       end
       
      % plot(mi.position(1),mi.position(2),'r.');
       
       ntargdist = sqrt((mi.position(1)-mi.target(1))^2+(mi.position(2)-mi.target(2))^2);
       
       if ntargdist < mi.targdist
          mi.targdist = ntargdist;
       elseif mi.age > .8 && ntargdist > mi.targdist %getting further away
           mi.hit = true;
           %delete(mi.graphic);
       end
               
       if ntargdist < .04
           mi.hit = true;
          
       end
       
      mi.trail = [mi.trail; mi.position, NaN];
      mi = plot(mi);  
   elseif mi.hit
       
       mi = explosion(mi);
       
       
   end
elseif mi.hit
   
    mi = explosion(mi);
end

end