function hit = buttonpress(abl,x,y)

if x <= 2 &&  x >= 1.95 && y <= 2 &&  y >= 1.95
    abl.launch = true;
  %  abl.button = fill(abl.bxcoords,abl.bycoords,'b');
    hit = true;
else
    hit = false;
end

end