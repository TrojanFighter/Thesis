function mi = cls_interceptor(px, py)
%Developed 08/2008
%Author David Robinson
%email cosmichero@gmail.com
mi.position = [px py];
mi.trail = [mi.position,NaN];
mi.velocity = [0 .03];
mi.heading  = 0;
mi.bodyaxis = [0 1;-1 0]; %start pointing up
mi.hit = false;
mi.broken = false;
mi.launched = false;
mi.target = [0 0];
mi.targdist = 10;
mi.age = 0;
mi.heading = 0;
mi.thrusting = false;
mi.explosion = [0,0,NaN]; %explosion time, radious, graphics handle
mi.stop = false;

fv.vertices = [-.025,  .025;...  %pt 1 border
               -.0125  .012;... %pt 2 border
               -.025   .012;... %pt 3
                .015   .012;... %pt 4 border
                .0375  .0   ;... %pt 5 border
                .015  -.012;... %pt 6 border
               -.0125 -.012;... %pt 7 border
               -.025  -.025; ... %pt 8 border
               -.025  -.012];   %pt 9
           
fv.faces = [1 2 3;... % left fin
            4 5 6;... % nose
            7 8 9;... %right fin
            3 4 9;... %body
            4 6 9];   %body
        
fv.FaceVertexCData = [0 0 0;
                      0 0 0;
                      0 0 0;
                      .5 .5 .5;
                      .5 .5 .5];
                  
mi.patchobj = fv;
rotmatrix = [cos(pi/2) sin(pi/2); -sin(pi/2) cos(pi/2)];
fv.vertices = [rotmatrix'*fv.vertices']';
fv.vertices(:,1) = fv.vertices(:,1) + mi.position(1);
fv.vertices(:,2) = fv.vertices(:,2) + mi.position(2);
mi.graphic = patch(fv);
set(mi.graphic,'EdgeColor',[.5 .5 .5])
set(mi.graphic,'FaceColor','flat');
% set(mi.graphic,'CData',fv.FaceVertexCData);
mi = class(mi,'cls_interceptor');

end