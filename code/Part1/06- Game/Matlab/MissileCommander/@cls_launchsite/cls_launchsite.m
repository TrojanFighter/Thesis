function ls  = cls_launchsite(x)
%Developed 08/2008
%Author David Robinson
%email cosmichero@gmail.com
global interceptors;

ls.xoffset = x;
ls.explosion = [0,0,NaN]; %explosion time, radious, graphics handle
ls.xcoords = [x x x+.2 x+.2];
ls.ycoords = [.6 .4 .4 .6];
ls.hitcount = 0;
ls.damagelines = [];
ls.launchfacility = fill(ls.xcoords,ls.ycoords,'y');%'FaceColor',[.5 .5 .5]);
set(ls.launchfacility,'FaceColor',[.5 .5 .5]);
ls.stop = false;
ls.hit = false;

rtheta = pi:.01:2*pi;
rx = cos(rtheta);
ry = sin(rtheta)+1;
rx = [rx,-1]*0.035+x;
ry = [ry,1]*0.035+.65;
ls.rpts = [rx',ry'];

ax = [-.1 -.1 .1 .1]*0.035+x;
ay = ([1 0 0 1]+1)*0.035+.65;
ls.apts = [ax',ay'];

sx = [-.125 -.125 .125 .125]*0.035+x;
sy = [0 -2 -2 0]*0.035+.65;


ls.atena = fill(ax,ay,'k');
ls.stand = fill(sx,sy,'k');
ls.radardish = fill(rx,ry,'r');



amunition = cell(1,1);
for i = 1:10
    offset = x+.02 +.015*i; 
    amunition{i} = cls_interceptor(offset,.6);
    ls.interceptors(i) = i;
end

interceptors = [interceptors,amunition];
ls.maxcount = max(size(interceptors));
ls.launchcount = max(size(interceptors))-9;



ls = class(ls,'cls_launchsite');


end