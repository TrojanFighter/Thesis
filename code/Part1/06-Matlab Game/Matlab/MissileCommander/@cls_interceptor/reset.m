function mi = reset(mi)

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

end