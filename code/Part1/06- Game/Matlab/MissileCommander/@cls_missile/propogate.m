function m = propogate(m)
global dt
global surface
global launchfacility1
global launchfacility2
global launchfacility3
global round
if m.launch && ~m.stop
   if ~m.hit %if not hit then continue flying along
      m.age = m.age+dt;
      m.position(1) = m.position(1)+ m.velocity(1)*dt;
      m.position(2) = m.position(2)+ m.velocity(2)*dt;
      m.velocity(2) = m.velocity(2) - .03*dt*((round+1)/2);
      m.trail = [m.trail; m.position, NaN];
      m = plot(m);
      if detectcollison(surface,m.position(1),m.position(2))||m.position(2) < 0 
          m.hit = true;
          
          surface = deformsurface(surface,m.position(1),m.position(2));
      end
      if detectcollison(launchfacility1,m.position(1),m.position(2))    
         m.hit = true;
         
         launchfacility1 = hit(launchfacility1);
      end
      if detectcollison(launchfacility2,m.position(1),m.position(2))    
         m.hit = true;
         
         launchfacility2 = hit(launchfacility2);
      end
      if detectcollison(launchfacility3,m.position(1),m.position(2))    
         m.hit = true;
         
         launchfacility3 = hit(launchfacility3);
      end
   elseif m.hit %if it is hit then explode
       
       m = explosion(m); 
%        if m.stop && detectcollison(surface,m.position(1),m.position(2))
%            surface = deformsurface(surface,m.position(1),m.position(2));
%        end
   end
end

% global dt
%   %while ~m.stop
%   if m.launch && ~m.stop 
%       
%       m.age = m.age+dt;
%       m.position(1) = m.position(1)+ m.velocity(1)*dt;
%       m.position(2) = m.position(2)+ m.velocity(2)*dt;
%       m.velocity(2) = m.velocity(2) - .1*dt;
%       m.trail = [m.trail; m.position, NaN];
%       m = plot(m);
%       
%       if m.position(2) < -1
%          m.stop = true;
%          delete(m); 
%       end
%       
%      
%           
%           
%       
%   end

end