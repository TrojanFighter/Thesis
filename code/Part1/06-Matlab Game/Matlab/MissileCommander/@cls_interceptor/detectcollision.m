function collision = detectcollision(mi,x,y)




if ~mi.hit 
    pts = get(mi.graphic,'vertices');
    if x <= max(pts(:,1)) && x >= min(pts(:,1)) &&...
         y <= max(pts(:,2)) && y >= min(pts(:,2))   
     
         collision = true;
    else
        collision = false;
    end
else
    distance = sqrt((mi.position(1)-x)^2+(mi.position(2)-y)^2);
    if distance <= mi.explosion(2)
        collision = true;
    else
        collision = false;
    end
end

if mi.launched == false
    collision = false;
end


end