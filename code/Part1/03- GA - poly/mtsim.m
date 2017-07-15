function f=mtsim(x)
global XNT
XNT=[0:20:100;x(1:6)]';

sim('MissileGuidanceMainBlockTrail1')
Ammax=max(XNC);
misD=Rtm(length(Rtm));
w1=10^5;
w2=1;
[f]=10^10/(w1*(misD).^2+w2*(Ammax).^2);