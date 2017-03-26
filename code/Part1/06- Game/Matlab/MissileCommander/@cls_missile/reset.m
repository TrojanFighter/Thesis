function m = reset(m)

delete(m);
m.age = 0;
m.explosion = [0,0,NaN]; %explosion time, radious, graphics handle
m.hit = false;
m.launch = false;
m.stop = false;
m.position = [rand()*2,2.5];
m.trail = [m.position,NaN];
m.shotdown = false;


end