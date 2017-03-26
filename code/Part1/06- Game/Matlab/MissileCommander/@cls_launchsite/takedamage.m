function ls = takedamage(ls)
global interceptors;
if  ls.hitcount > 0 && ls.hitcount < 3 && ls.hit
    xpts = rand(5,1)*.2 + ls.xoffset;
    xpts = sort(xpts);
    ypts = rand(5,1)*.2 + .4;
    ls.damagelines(ls.hitcount) = line(xpts,ypts,'Color','k');
    ls.hit = false;
elseif ls.hitcount >= 3  && ~ls.stop
    for i = ls.maxcount-9:ls.maxcount
       interceptors{i}=hit(interceptors{i}); 
    end    
    ls = explosion(ls);
end