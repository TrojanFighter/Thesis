function m = explosion(m)
    global dt
    global score
    global explosionplayer;
    m.explosion(1) = m.explosion(1) + dt;
    if m.explosion(1) == dt
       % wavplay(explosionfx{1},explosionfx{2},'async')
       delete(m)
       stop(explosionplayer);
       play(explosionplayer);
    end
    m.explosion(2) = m.explosion(2).^2 + dt;
    if ~isnan(m.explosion(3))
        delete(m.explosion(3));
    end
    m.explosion(3) = plotcircle(m.explosion(2),m.position(1),m.position(2));
    if 1-m.explosion(1) > 0
        set(m.explosion(3),'facealpha',1-m.explosion(1));
        set(m.explosion(3),'edgecolor','y');
    else
        delete(m.explosion(3));
        m.position = [0 0];
        m.stop = true;
        if m.shotdown == true
            score = score + 100;
        end
    end
end


function h = plotcircle(R,xd,yd)
t = 0:.1:2*pi;
y = R*cos(t)+yd;
x = R*sin(t)+xd;
hold on;
h = fill(x,y,[1,rand()*.5,rand()*.5]);
end %function