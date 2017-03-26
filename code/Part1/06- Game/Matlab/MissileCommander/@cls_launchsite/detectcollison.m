function inside = detectcollison(ls,x,y)
xp = ls.xcoords;
yp = ls.ycoords;
edgecount = 0;

%http://geometryalgorithms.com/Archive/algorithm_0103/algorithm_0103.htm

%1. an upward edge includes its starting endpoint, and excludes its final endpoint; 
% 
%2. a downward edge excludes its starting endpoint, and includes its final endpoint; 
% 
%3. horizontal edges are excluded; and 
% 
%4. the edge-ray intersection point must be strictly right of the point P


if x > min(xp) && x < max(xp) && y > min(yp) && y < max(yp)
    for i = 1:max(size(yp))
        if i == max(size(yp))
            edge = [xp(i) xp(1) yp(i) yp(1)];
        else
            edge = [xp(i) xp(i+1) yp(i) yp(i+1)];
        end
        if x < max(edge(1:2)) %testing for crossings to the 'right' of x
            %test the kind of edge
            if edge(3)== edge(4)
                type = 'horizontal';
            elseif edge(3) < edge(4)
                type = 'upward';
            elseif edge(3) > edge(4)
                type = 'downward';
            end

            switch type
                case 'horizontal'
                    %does not count
                case 'upward'
                    if y >= edge(3) && y < edge(4) %include begining point
                        edgecount = edgecount+1;
                    end
                case 'downward'
                    if y < edge(3) && y >= edge(4) %include end point
                        edgecount = edgecount+1;
                    end
            end
        end   
    end
    if mod(edgecount,2)==0
        inside = false;
    else
        inside = true;
    end
else
    inside = false; 
end
    



    
end



    



    
