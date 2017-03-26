function ls = explosion(ls)
    global dt
    global explosionplayer;
    ls.explosion(1) = ls.explosion(1) + dt;
    if ls.explosion(1) == dt
       % wavplay(explosionfx{1},explosionfx{2},'async')
       stop(explosionplayer);
       play(explosionplayer);
    end
    ls.explosion(2) = ls.explosion(2).^2 + dt*2;
    if ~isnan(ls.explosion(3))
        delete(ls.explosion(3));
    end
    ls.explosion(3) = plotcircle(ls.explosion(2),ls.xoffset+.1,.4);
    if 1-ls.explosion(1) > 0
        set(ls.explosion(3),'facealpha',1-ls.explosion(1)/2);
        set(ls.explosion(3),'edgecolor','r');
    else
        delete(ls.explosion(3));
        ls.stop = true;
        for i = 1:max(size(ls.damagelines))
        delete(ls.damagelines(i));
        end
        delete(ls.launchfacility);
        delete(ls.atena); 
        delete(ls.stand); 
        delete(ls.radardish); 
        xtop = [0:.04:.2]+ls.xoffset;
        ytop = rand(size(xtop))*.1 + .4;
        ls.xcoords = [xtop,ls.xoffset+.2,ls.xoffset];
        ls.ycoords = [ytop,.4,.4];
        ls.launchfacility = fill(ls.xcoords,ls.ycoords,'y');%'FaceColor',[.5 .5 .5]);
        set(ls.launchfacility,'FaceColor',[.5 .5 .5]);  
    end
end


function h = plotcircle(R,xd,yd)
t = 0:.1:2*pi;
y = R*cos(t)+yd;
x = R*sin(t)+xd;
hold on;
h = fill(x,y,[0,0,1]);
end %function