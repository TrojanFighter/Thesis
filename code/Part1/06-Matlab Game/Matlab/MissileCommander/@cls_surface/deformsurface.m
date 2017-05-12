function s = deformsurface(s,x,y)
xp = s.x;
yp = s.t;



li = find(xp < x-.05);
gi = find(xp > x+.05);

 

if ~isempty(li) && ~isempty(gi)
    
    gi = gi(1);
    
    li = li(end);
    
    
    xstart = xp(li)-x;
    xend = xp(gi)-x;
    ystart = yp(li)-y;
    yend = yp(gi)-y;
    
    
    if ystart >= 0
        startangle = acos(xstart/sqrt(ystart^2+xstart^2));
    else
        startangle = 2*pi - acos(xstart/sqrt(ystart^2+xstart^2));
    end
    
    if yend >= 0
        endangle = acos(xend/sqrt(yend^2+xend^2));
    else
        endangle = 2*pi - acos(xend/sqrt(yend^2+xend^2));
    end
    
    if endangle < startangle 
        endangle = endangle + 2*pi;
    end
    
    theta = startangle:.1:endangle;
    if startangle == 0 && endangle == 0
        theta = 0:.1:pi;
    end
    
    xc = .05*cos(theta)+x;
    yc = .05*sin(theta)+y;
    
    
    
    xp = [xp(1:li),xc,xp(gi:end)];
    yp = [yp(1:li),yc,yp(gi:end)];
    
    delete(s.f);        
    s.f = fill([-.5,xp,2.5],[min(yp)-1,yp,min(yp)-1],'g');
   
    
    
    
    
    
%     distance = sqrt((xp(i)-x)^2 + (yp(i)-y)^2);
%     if distance < .2
%         if distance ~= 0
%             if .002/distance > max
%                 deform = max;
%             else
%                 deform = .002/distance;
%             end
%         elseif distance == 0
%             deform = max
%         end
%         
%         if xp(i)-x >= 0
%             xp(i) = xp(i)+ deform;
%         elseif xp(i)-x < 0
%             xp(i) = xp(i)- deform;
%         end        
%         yp(i) = yp(i) - deform*2;
%         delete(s.f);
%         s.f = fill([0,xp,2],[min(yp)-1,yp,min(yp)-1],'g');    
%     end



s.x = xp;
s.t = yp;
end

end




