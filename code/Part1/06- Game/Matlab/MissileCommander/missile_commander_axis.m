function missile_commander_axis
%Developed 08/2008
%Author David Robinson
%email cosmichero@gmail.com

%set up a figure to play the misssile commander Safegaurd game


h = figure;
set(h,'Position',[50 50 840 630]);
a = axes;
axis([0 2 0 2]);
set(a,'Color','c');
set(a,'XTick',[],'YTick',[]);
set(a,'NextPlot','add');
title('Safeguard: Round One','Color','b','FontSize',18)

%call these functions when these figure events occur
set(h,'WindowButtonMotionFcn',@RotateRadar)
set(h,'WindowButtonDownFcn',@setheading);
%set(h,'KeyPressFcn',@keypressinfo) %buggy


end %function





function RotateRadar(src,evnt) %makes the radars point at the mouse
 global launchfacility1;
 global launchfacility2;
 global launchfacility3;
 
 c = get(src,'CurrentAxes');
 m = get(c,'CurrentPoint');
launchfacility1 = AimRadar(launchfacility1,m(1,1),m(1,2));
launchfacility2 = AimRadar(launchfacility2,m(1,1),m(1,2));
launchfacility3 = AimRadar(launchfacility3,m(1,1),m(1,2));
end

function setheading(src,evnt)

global launchfacility1;
global launchfacility2;
global launchfacility3;
global activelauncher;
global airforce1;
global specialused;
c = get(src,'CurrentAxes');
pt = get(c,'CurrentPoint');

%first check to see if any of the clickable objects were clicked
% (launch facilities or the special weapon button)
if detectcollison(launchfacility1,pt(1,1),pt(1,2))
    activelauncher = 1;
    highlight(launchfacility1,'on');
    highlight(launchfacility2,'off');
    highlight(launchfacility3,'off');
elseif detectcollison(launchfacility2,pt(1,1),pt(1,2))
    activelauncher = 2;
    highlight(launchfacility1,'off');
    highlight(launchfacility2,'on');
    highlight(launchfacility3,'off');
elseif detectcollison(launchfacility3,pt(1,1),pt(1,2))
    activelauncher = 3;
    highlight(launchfacility1,'off');
    highlight(launchfacility2,'off');
    highlight(launchfacility3,'on');
elseif  buttonpress(airforce1,pt(1,1),pt(1,2))
    if specialused == false;
        specialused = true;
        airforce1 = takeoff(airforce1);
    end
else %launch a missile from the active launcher
    switch activelauncher 
        case 1
           launchfacility1 = launchinterceptor(launchfacility1,pt(1,1),pt(1,2)); 
       case 2
           launchfacility2 = launchinterceptor(launchfacility2,pt(1,1),pt(1,2)); 
       case 3
           launchfacility3 = launchinterceptor(launchfacility3,pt(1,1),pt(1,2));
    end
end


end

function keypressinfo(src1,evnt1)
global activelauncher
global launchfacility1;
global launchfacility2;
global launchfacility3;

switch evnt1.Key
    case '1'
       activelauncher = 1;
       highlight(launchfacility1,'on');
       highlight(launchfacility2,'off');
       highlight(launchfacility3,'off');
    case '2'
        activelauncher = 2;
        highlight(launchfacility1,'off');
        highlight(launchfacility2,'on');
        highlight(launchfacility3,'off');
    case '3'
        activelauncher = 3;
        highlight(launchfacility1,'off');
        highlight(launchfacility2,'off');
        highlight(launchfacility3,'on');
end
        

end



