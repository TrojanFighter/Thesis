function mi = explosion(mi)
    global dt
    global explosionplayer;
    
    
    if ~mi.stop
    mi.explosion(1) = mi.explosion(1) + dt;
    if mi.explosion(1) == dt
       % wavplay(explosionfx{1},explosionfx{2},'async')
       delete(mi.graphic)
       stop(explosionplayer);
       play(explosionplayer);
    end
    mi.explosion(2) = mi.explosion(2).^2 + dt*1.5;
    if ~isnan(mi.explosion(3)) 
       delete(mi.explosion(3)); 
    end
    mi.explosion(3) = plotcircle(mi.explosion(2),mi.position(1),mi.position(2));
    if 1-mi.explosion(1) > 0
        set(mi.explosion(3),'facealpha',1-mi.explosion(1)/5);
        set(mi.explosion(3),'edgecolor','r');
    else
        mi.stop = true;
        mi.position = [0,0];
        delete(mi.explosion(3));
        delete(mi);
    end
    end
end


function h = plotcircle(R,xd,yd)
t = 0:.1:2*pi;
y = R*cos(t)+yd;
x = R*sin(t)+xd;
hold on;
h = fill(x,y,[1,1,rand()*.5]);
end %function