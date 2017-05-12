function m = cls_missile()
%Developed 08/2008
%Author David Robinson
%email cosmichero@gmail.com



m.age = 0;
m.explosion = [0,0,NaN]; %explosion time, radious, graphics handle
m.hit = false;
m.launch = false;
m.stop = false;
m.position = [rand()*2,2.5];
m.trail = [m.position,NaN];
m.shotdown = false;
if m.position(1) >= 1
    dir = -1;
else
    dir = 1;
end
    
m.velocity = [dir*rand()*.3,-rand()*.2];

fv.vertices = [-.0125 .0125;...
                .0100 .0100;...
                .0375  0;  ...
                .0100 -.0100;...
               -.0125 -.0125];
           
fv.faces = [1 2 5;... %body
            2 4 5;... %body
            2 3 4];... %nose
        
fv.FaceVertexCData = [0 0 0;...
                      0 0 0;...
                      1 0 0];
                      
                      
                  
m.patchobj = fv;

m.patchobj = fv;
rotmatrix = [cos(pi/2) sin(pi/2); -sin(pi/2) cos(pi/2)];
fv.vertices = [rotmatrix'*fv.vertices']';
fv.vertices(:,1) = fv.vertices(:,1) + m.position(1);
fv.vertices(:,2) = fv.vertices(:,2) + m.position(2);
m.graphic = patch(fv);
set(m.graphic,'EdgeColor',[0 0 0])
set(m.graphic,'FaceColor','flat');





m = class(m,'cls_missile');





end