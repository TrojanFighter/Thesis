function abl = takeoff(abl)
 if abl.launch == false
  abl.launch = true;
  abl.button = fill(abl.bxcoords,abl.bycoords,'b');
 end
end