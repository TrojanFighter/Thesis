function [RTM,XNCmax] = PNeq (polyORtrapORtrapSymm,polydegree, c , draw)

%TWO-DIMENSIONAL TACTICAL MISSILE-TARGET ENGAGEMENT SIMULATION
%solving using the second-order Runge–Kutta numerical integration technique


%% Simulation inputs
G=32.2;
n=0; %counter on points
VM = 3000.; %magnitude of the missile velocity [ft/sec]
VT = 1000.; %magnitude of the target velocity [ft/sec]
HEDEG = 0.; %heading error (in degrees)       [second error source]
XNP = 5.; %effective navigation ratio [dimensionless]
RM1 = 0.0; %initial location of the missile (in the dowenrange axis)[ft]
RM2 = 1.0; %initial location of the missile (in the crossrange axis)[ft]
RT1 = 12000.; %initial location of the target (in the dowenrange axis)[ft]
RT2 = 10000.; %initial location of the target (in the crossrange axis)[ft]
BETA=0.; %angle between the target velocity vector and the doewrange axis

%% Differential equations

VT1=-VT*cos(BETA); %component of the target velocity vector in downrange axis
VT2=VT*sin(BETA); %component of the target velocity vector in crossrange axis
HE=HEDEG/57.3; %heading error (in rad)

T=0.; %time
S=0.;
RTM1=RT1-RM1; %relative missile-target separations (downrange)
RTM2=RT2-RM2; %relative missile-target separations (crossrange)
RTM=sqrt(RTM1*RTM1+RTM2*RTM2); %relative separation between missile and target
XLAM=atan2(RTM2,RTM1); %line-of-sight angle
XLEAD=asin(VT*sin(BETA+XLAM)/VM); %missile lead angle
THET=XLAM+XLEAD;
VM1=VM*cos(THET+HE); %missile velocity component (in downrange)
VM2=VM*sin(THET+HE); %missile velocity component (in crossrange)
VTM1 = VT1 - VM1; %relative velocity component (in downrange)
VTM2 = VT2 - VM2; %relative velocity component (in crossrange)
VC=-(RTM1*VTM1 + RTM2*VTM2)/RTM %closing velocity

iter=0;

while VC >= 0 %terminate the programme when the velocity chnges its sign
    %means that the separation between the missile and target is a minimum
    iter=iter+1;
    if RTM < 1000
        H=.0002; %step size made smaller near the end of the flight
        %(to accurately capture the magnitude of the miss distance)
    else
        H=.01; %step of the most of the flight (except the end)
    end
    
    if strcmp(polyORtrapORtrapSymm,'poly')
        T1 = polydegree-1 : -1 :0;
        XNT=sum(c.*(T.^T1)); %target acceleration (target manuver)[first error source]
    elseif strcmp(polyORtrapORtrapSymm,'trap')
        XNT=trapacc(c(1),c(2),c(3),c(4),c(5),T);
    elseif strcmp(polyORtrapORtrapSymm,'trapSymm')
        XNT=trapacc(c(1),c(2),c(3),c(4),c(5),T);
    end
    
    
    %---------------
    
    BETAOLD=BETA;
    RT1OLD=RT1;
    RT2OLD=RT2;
    RM1OLD=RM1;
    RM2OLD=RM2;
    VM1OLD=VM1;
    VM2OLD=VM2;
    STEP=1;
    FLAG=0;
    while STEP <=1
        if FLAG==1
            STEP=2;
            BETA=BETA+H*BETAD;
            RT1=RT1+H*VT1;
            RT2=RT2+H*VT2;
            RM1=RM1+H*VM1;
            RM2=RM2+H*VM2;
            VM1=VM1+H*AM1;
            VM2=VM2+H*AM2;
            T=T+H;
        end
        RTM1=RT1-RM1;
        RTM2=RT2-RM2;
        RTM=sqrt(RTM1*RTM1+RTM2*RTM2);
        VTM1=VT1-VM1;
        VTM2=VT2-VM2;
        VC=-(RTM1*VTM1+RTM2*VTM2)/RTM;
        XLAM=atan2(RTM2,RTM1);
        XLAMD=(RTM1*VTM2-RTM2*VTM1)/(RTM*RTM);
        XNC=XNP*VC*XLAMD; %command acceleration of the rocket[proportional navigation guidance law]
        if XNC>60*G
            XNC=60*G;
        end
        xnc_max(iter)=XNC;
        AM1=-XNC*sin(XLAM);
        AM2=XNC*cos(XLAM);
        VT1=-VT*cos(BETA);
        VT2=VT*sin(BETA);
        BETAD=XNT/VT;
        FLAG=1;
        
        
    end
    FLAG=0;
    BETA=.5*(BETAOLD+BETA+H*BETAD);
    RT1=.5*(RT1OLD+RT1+H*VT1);
    RT2=.5*(RT2OLD+RT2+H*VT2);
    RM1=.5*(RM1OLD+RM1+H*VM1);
    RM2=.5*(RM2OLD+RM2+H*VM2);
    VM1=.5*(VM1OLD+VM1+H*AM1);
    VM2=.5*(VM2OLD+VM2+H*AM2);
    S=S+H;
    if S>=.09999
        S=0.;
        n=n+1;
        ArrayT(n)=T;
        ArrayXNTG(n)=XNT/32.2;
        ArrayRT1(n)=RT1;
        ArrayRT2(n)=RT2;
        ArrayRM1(n)=RM1;
        ArrayRM2(n)=RM2;
        ArrayXNCG(n)=XNC/32.2;
        ArrayRTM(n)=RTM;
        
        
        
    end
    
end

RTM; %miss distance
XNCmax=max(xnc_max);

%%----outputs----
if draw==1 %to viwe all figurs (runs)
    h = figure('position',[100, 100, 1000, 750]); set(gcf,'color','w'); set(gca,'FontSize',24);
    plot(ArrayRT1,ArrayRT2,':',ArrayRM1,ArrayRM2,'*r','LineWidth',3),grid on
    h_legend = legend('Target trajectory','Missile trajectory','Location','NorthWest');
    set(h_legend,'FontSize',14); set(gca,'FontSize',14);
    title('Two-dimensional tactical missile-target engagement simulation')
    xlabel('Downrange (Ft) ','FontSize', 16)
    ylabel('Altitude or crossrange (Ft)','FontSize', 16)
    saveTightFigure(h,'trajectoryTRAPSymm.pdf')
    
    hold on   %-------------------------------
    
    h = figure('position',[100, 100, 1000, 750]); set(gcf,'color','w'); set(gca,'FontSize',24);
    plot(ArrayT,ArrayXNCG,'LineWidth',3),grid on
    %title('Two-dimensional tactical missile-target engagement simulation')
    xlabel('Time (sec)','FontSize', 16)
    ylabel('Acceleration of missle (G)','FontSize', 16)
    saveTightFigure(h,'MissileAccelerationTRAPSymm.pdf')
    
    hold on
    
    h = figure('position',[100, 100, 1000, 750]); set(gcf,'color','w'); set(gca,'FontSize',24);
    plot(ArrayT,ArrayXNTG,'LineWidth',3),grid on
    %title('Two-dimensional tactical missile-target engagement simulation')
    xlabel('Time (sec)','FontSize', 16)
    ylabel('Acceleration of Target (G)','FontSize', 16)
    saveTightFigure(h,'TargetAccelerationTRAPSymm.pdf')
end