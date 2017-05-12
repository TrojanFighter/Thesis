function  s = cls_surface()
%Developed 08/2008
%Author David Robinson
%email cosmichero@gmail.com
s.x = -.5:.01:2.5;
s.x = round(s.x*100)/100;
y = sin(s.x*8)*.125+.5;
z = rand(size(y))*.125-.0625;
s.t = y+z;
s.p = polyfit(s.x,s.t,12);
s.t = polyval(s.p,s.x);

%make a spot for 3 launch facilities



side1y = .6:-.01:.4;
side1x = ones(size(side1y));
bottomx = 0:.01:.2;
bottomy = ones(size(bottomx));


i1 = find(s.x == .1);
i2 = find(s.x == .3);
s.x = [s.x(1:i1),[side1x*.1,bottomx+.1,side1x*.3],s.x(i2:end)];
s.t = [s.t(1:i1-1),[.6,side1y,bottomy*.4,side1y,.6],s.t(i2+1:end)];

i1 = find(s.x == .9);
i2 = find(s.x == 1.1);
s.x = [s.x(1:i1),[side1x*.9,bottomx+.9,side1x*1.1],s.x(i2:end)];
s.t = [s.t(1:i1-1),[.6,side1y,bottomy*.4,side1y,.6],s.t(i2+1:end)];


i1 = find(s.x == 1.7);
i2 = find(s.x == 1.9);
s.x = [s.x(1:i1),[side1x*1.7,bottomx+1.7,side1x*1.9],s.x(i2:end)];
s.t = [s.t(1:i1-1),[.6,side1y,bottomy*.4,side1y,.6],s.t(i2+1:end)];



s.f = fill([-.5,s.x,2.5],[min(s.t)-1,s.t,min(s.t)-1],'g');
s = class(s,'cls_surface');
end