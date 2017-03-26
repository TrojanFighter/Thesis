function abl = propogate(abl)
global missiles
global phaserplayer
global jetplayer;
global hl2dachfplayer;
if abl.launch == true
if abl.position(1) <= 2.5
if abl.position(1) == -.5
    play(jetplayer)
    play(hl2dachfplayer)
end
delete(abl);
abl.position(1) = abl.position(1) + .05;
abl = plotABL(abl);

else
   abl.launch = false;
   abl.position(1) = -.5;
   abl.button = fill(abl.bxcoords,abl.bycoords,'k');
end

if abl.lasertime == 0
    for i = 1:max(size(missiles))
        mpos = get(missiles{i},'Position');
        if mpos(1) > abl.position(1)
        distance = sqrt((mpos(1)-abl.nose(1))^2 + (mpos(2)-abl.nose(2))^2);
        if distance <= .8 && get(missiles{i},'hit')==false
            abl.targetindex = i;
            abl.lasertime = 1;
            break
        end
        end
    end
elseif abl.lasertime >= 1 &&  abl.lasertime < 4
    if abl.lasertime == 1
        play(phaserplayer);
    end
    try 
        delete(abl.laser);
    end
    
    mpos = get(missiles{abl.targetindex},'Position');
    abl.nose(1);
    abl.nose(2);
    abl.laser = line([abl.nose(1),mpos(1)],[abl.nose(2),mpos(2)],'Color','r');
    abl.lasertime = abl.lasertime + 1;
    missiles{abl.targetindex} = lasered(missiles{abl.targetindex});
elseif abl.lasertime == 4 
    missiles{abl.targetindex}=hit(missiles{abl.targetindex});
    try 
        delete(abl.laser)
    end
    abl.lasertime = abl.lasertime + 1;
elseif abl.lasertime >= 5 && abl.lasertime <= 8
    abl.lasertime = abl.lasertime + 1;
    if abl.lasertime == 7
        abl.lasertime = 0;
    end   
end
        


end
end
