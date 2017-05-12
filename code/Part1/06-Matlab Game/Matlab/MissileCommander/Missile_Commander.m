function Missile_Commander
%Developed 08/2008
%Author David Robinson
%email cosmichero@gmail.com

%declare global variables
global score
global explosionplayer;
global hl2dachfplayer;
global phaserplayer;
global jetplayer;
global launchplayer;
global dt;                %Delta time -time step used throughout the game
global surface;
global launchfacility1;
global launchfacility2;
global launchfacility3;
global interceptors;
global activelauncher;
global missiles;
global airforce1;
global round;

global specialused;




%call function to set up the axis
missile_commander_axis

%initialize globals:
specialused = false;
score = 0;
activelauncher = 1;
interceptors = [];
[y, fs] = audioread('explosion-01.wav');
explosionplayer = audioplayer(y, fs);
[y, fs] = audioread('missile.wav');
launchplayer = audioplayer(y, fs);
[y, fs] = audioread('shipphaser1.wav');
phaserplayer = audioplayer(y, fs);
[y, fs,] = audioread('jetflyby3.wav');
jetplayer = audioplayer(y, fs);
[y, fs] = audioread('hl2dachf.wav');
hl2dachfplayer = audioplayer(y, fs);
dt = .1;
surface = cls_surface;
launchfacility1  = cls_launchsite(.1);
highlight(launchfacility1,'on');
launchfacility2  = cls_launchsite(.9);
launchfacility3  = cls_launchsite(1.7);
airforce1 = cls_abl();
missiles = cell(1,1);
for i = 1:10
  missiles{i} = cls_missile;
end

launchcount = 0; %used to determine which missile to launch next

for round = 1:3 %game lasts for 3 rounds
    
    %Update the title with the current round
    if round == 1
    title('Safeguard: Round One','Color','b','FontSize',18)
    elseif round == 2
        title('Safeguard: Round Two','Color','b','FontSize',18)
    elseif round == 3
        title('Safeguard: Round Three','Color','b','FontSize',18)
    end

    for i = 1:dt:65 %go to 65 seconds per round

         if i < 55 %don't launch a missile if there are only 10 seconds left in the round
             if mod(i,1)== 0 && launchcount <=10 %every time i is an integer launch a missile
                launchcount = launchcount + 1;
                missiles{launchcount}= launchmissile(missiles{launchcount});
                if launchcount == 10
                    launchcount = 0;
                end
             end
         end %if i < 60



        for j=1:10 %propogate all the missiles
           missiles{j}=propogate(missiles{j}); 
        end

        for j=1:max(size(interceptors))%propogate all the interceptors
           interceptors{j}=propogate(interceptors{j});
        end

        %determine if any interceptors hit any of the missiles
        for k = 1:max(size(interceptors))
            for j = 1:10
                   position = get(missiles{j},'Position');
                   if detectcollision(interceptors{k},position(1),position(2))
                      missiles{j}=hit(missiles{j});
                      missiles{j}=shotdown(missiles{j});
                      interceptors{k}=hit(interceptors{k});
                   end     
            end
        end

        %update see if any of the launch facilities take damage and propogate
        %airforce1
        launchfacility1 = takedamage(launchfacility1);
        launchfacility2 = takedamage(launchfacility2);
        launchfacility3 = takedamage(launchfacility3);
        airforce1 = propogate(airforce1);

        %check to see if any of the missiles have stopped if so reintialize
        %them
        for j=1:10
           if get( missiles{j},'stopped')
                missiles{j}= cls_missile; 
           end
        end

        %update the score
        xlabel(['Score: ',num2str(score)],'Color','r','FontSize',15);
        pause(dt);

        %if all the launch facilites are destroyed stop
        if isdestroyed(launchfacility1) && isdestroyed(launchfacility2) && isdestroyed(launchfacility3)
           break
        end


    end %for i = 1:dt:65 

%reset for next round
specialused = false;
airforce1=changebuttoncolor(airforce1,'r');

%delete the missiles
for j=1:10
   delete(missiles{j}); 
end

%get rid of the interceptors
for j=1:max(size(interceptors))
   delete(interceptors{j});
end

%reinitialize the missiles
for i = 1:10
  missiles{i} = cls_missile;
end

%reset the interceptors and reload facilities
interceptors = [];
launchfacility1 = reload(launchfacility1);
launchfacility2 = reload(launchfacility2);
launchfacility3 = reload(launchfacility3);

%if the all the facilities are destroyed end the game
if isdestroyed(launchfacility1) && isdestroyed(launchfacility2) && isdestroyed(launchfacility3)
   break
end

end %round


%round loop completes when you either survive or all the facilities are
%destroyed
if isdestroyed(launchfacility1) && isdestroyed(launchfacility2) && isdestroyed(launchfacility3)
   title('Safeguard: Game Over All hope is lost','Color','r','FontSize',18);
else
   title('Safeguard: Congrats You Survived','Color','b','FontSize',18);
end

 
end

