function f=mtsim_trap(d) %these are the values of the time intervals between the times to change maneuvers
global XNT
G=32.2;
max_acc = 10* G;

d = abs(d);
% if d==0
%    d=d+0.001 
% end
t(1) = 0;
for i = 2:length(d)
    t(i) = t(i-1) + d(i);
end
t(length(d)+1) = 60;

XNT=[t;[0 max_acc max_acc 0 -max_acc -max_acc 0]]'; %here the t is the time values that we want to optimize

sim('MissileGuidanceMainBlock22')
Ammax=max(XNC);
misD=Rtm(length(Rtm));
w1=10^5;
w2=1;
[f]=10^10/(w1*(misD).^2+w2*(Ammax).^2);