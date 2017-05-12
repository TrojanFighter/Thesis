function ls = AimRadar(ls,x,y)

% ls.atena = fill(ax,ay,'b');
% ls.stand = fill(sx,sy,'y');
% ls.radardish = fill(rx,ry,'r');

%rotate the dish
if ls.hitcount < 1

pts = ls.rpts;

vx = x-ls.xoffset;
vy = y-.65;
theta = -acos(vx/sqrt(vx^2+vy^2))+pi/2;

 if vy < 0
      theta = 2*pi-theta; 
 end
if vy > 0
rx = [pts(:,1)-ls.xoffset]';
ry = [pts(:,2)-.65]';

pts = [cos(theta) sin(theta); -sin(theta) cos(theta)]*[rx;ry];
pts(1,:) = pts(1,:)+ls.xoffset;
pts(2,:) = pts(2,:)+.65;

delete(ls.radardish)
ls.radardish = fill(pts(1,:),pts(2,:),'r');
end

pts = ls.apts;
if vy > 0
rx = [pts(:,1)-ls.xoffset]';
ry = [pts(:,2)-.65]';

pts = [cos(theta) sin(theta); -sin(theta) cos(theta)]*[rx;ry];
pts(1,:) = pts(1,:)+ls.xoffset;
pts(2,:) = pts(2,:)+.65;

delete(ls.atena)
ls.atena = fill(pts(1,:),pts(2,:),'b');
end

end

end