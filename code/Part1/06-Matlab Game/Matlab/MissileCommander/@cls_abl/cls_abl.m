function abl = cls_abl()
%Developed 08/2008
%Author David Robinson
%email cosmichero@gmail.com
abl.position = [-.5 1];
abl.launch = false;
abl.graphicshandles = [];
abl.targetindex = 0;
abl.laser = [];
abl.nose = [];
abl.lasertime = 0;

abl.bxcoords = [1.95 1.95 2 2];
abl.bycoords = [2 1.95 1.95  2];
abl.button = fill(abl.bxcoords,abl.bycoords,'r');
abl.star = plot(1.975,1.975,'wp','markersize',10);


abl = class(abl,'cls_abl');
abl = plotABL(abl);
end