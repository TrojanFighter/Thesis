function varargout = game(varargin)
% GAME MATLAB code for game.fig
%      GAME, by itself, creates a new GAME or raises the existing
%      singleton*.
%
%      H = GAME returns the handle to a new GAME or the handle to
%      the existing singleton*.
%
%      GAME('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GAME.M with the given input arguments.
%
%      GAME('Property','Value',...) creates a new GAME or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before game_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to game_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help game

% Last Modified by GUIDE v2.5 15-Mar-2017 12:37:32

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @game_OpeningFcn, ...
    'gui_OutputFcn',  @game_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before game is made visible.
function game_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to game (see VARARGIN)

% Choose default command line output for game
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes game wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = game_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function Rt1_Callback(hObject, eventdata, handles)
% hObject    handle to Rt1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Rt1 as text
%        str2double(get(hObject,'String')) returns contents of Rt1 as a double


% --- Executes during object creation, after setting all properties.
function Rt1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Rt1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Rt2_Callback(hObject, eventdata, handles)
% hObject    handle to Rt2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Rt2 as text
%        str2double(get(hObject,'String')) returns contents of Rt2 as a double


% --- Executes during object creation, after setting all properties.
function Rt2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Rt2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function VT_Callback(hObject, eventdata, handles)
% hObject    handle to VT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of VT as text
%        str2double(get(hObject,'String')) returns contents of VT as a double


% --- Executes during object creation, after setting all properties.
function VT_CreateFcn(hObject, eventdata, handles)
% hObject    handle to VT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Rm1_Callback(hObject, eventdata, handles)
% hObject    handle to Rm1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Rm1 as text
%        str2double(get(hObject,'String')) returns contents of Rm1 as a double


% --- Executes during object creation, after setting all properties.
function Rm1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Rm1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Rm2_Callback(hObject, eventdata, handles)
% hObject    handle to Rm2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Rm2 as text
%        str2double(get(hObject,'String')) returns contents of Rm2 as a double


% --- Executes during object creation, after setting all properties.
function Rm2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Rm2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function VM_Callback(hObject, eventdata, handles)
% hObject    handle to VM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of VM as text
%        str2double(get(hObject,'String')) returns contents of VM as a double


% --- Executes during object creation, after setting all properties.
function VM_CreateFcn(hObject, eventdata, handles)
% hObject    handle to VM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function HEDEG_Callback(hObject, eventdata, handles)
% hObject    handle to HEDEG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of HEDEG as text
%        str2double(get(hObject,'String')) returns contents of HEDEG as a double


% --- Executes during object creation, after setting all properties.
function HEDEG_CreateFcn(hObject, eventdata, handles)
% hObject    handle to HEDEG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function XNP_Callback(hObject, eventdata, handles)
% hObject    handle to XNP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of XNP as text
%        str2double(get(hObject,'String')) returns contents of XNP as a double


% --- Executes during object creation, after setting all properties.
function XNP_CreateFcn(hObject, eventdata, handles)
% hObject    handle to XNP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%TWO-DIMENSIONAL TACTICAL MISSILE-TARGET ENGAGEMENT SIMULATION
%solving using the second-order Runge–Kutta numerical integration technique

%% Simulation inputs
% t = datestr(datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z'));
% push=str2num(get(handles.push,'string'));
% push=push+1;
% set(handles.push,'string',push)
% set(handles.time,'string',t)
dx=str2num(get(handles.dx,'string'));
dy=str2num(get(handles.dy,'string'));

RT1=str2num(get(handles.Rt1,'string'));
RT2=str2num(get(handles.Rt2,'string'))+dy;
RM1=str2num(get(handles.Rm1,'string'));
RM2=str2num(get(handles.Rm2,'string'));
VM=str2num(get(handles.VM,'string'));
VT=str2num(get(handles.VT,'string'));
XNP=str2num(get(handles.XNP,'string'));
XNT=str2num(get(handles.XNTT,'string'));
HEDEG=str2num(get(handles.HEDEG,'string'));

set(handles.Rt2,'string',RT2);
tic
n=0; %counter on points
%VM = 3000.; %magnitude of the missile velocity
%VT = 1000.; %magintude of the target velocity
%XNTT = 96.6; %target accelration (target manuver)[first error source]
%HEDEG = 0.0; %heading error (in degrees)       [second error source]
%XNP = 5.; %effictive navigation ratio
%RM1 = 0.; %initial location of the missile (in the dowenrange axis)
%RM2 = 9000.; %initial location of the missile (in the crossrange axis)
%RT1 = 12000.; %initial location of the target (in the dowenrange axis)
%RT2 = 10000.; %initial location of the target (in the crossrange axis)
BETA=0.; %angle between the target velocity vector and the doewrange axis

%% Diffrential equations

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
VC=-(RTM1*VTM1 + RTM2*VTM2)/RTM; %closing velocity

while VC >= 0 %terminate the programme when the velocity chnges its sign
    %means that the separation between the missile and target is a minimum
    if RTM < 1000
        H=.0002; %step size made smaller near the end of the flight
        %(to accurately capture the magnitude of the miss distance)
    else
        H=.01; %step of the most of the flight (except the end)
    end
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
        XNC=XNP*VC*XLAMD;
        XNC1=(XNP+20000000*XLAMD^2/VC)*VC*XLAMD;
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
        ArrayRT1(n)=RT1;
        ArrayRT2(n)=RT2;
        ArrayRM1(n)=RM1;
        ArrayRM2(n)=RM2;
        ArrayXNCG(n)=XNC/32.2;
        ArrayRTM(n)=RTM;
    end
end
RTM %miss distance

x=floor(VM/VT);

set(handles.Rm1,'string',ArrayRM1(x));
set(handles.Rm2,'string',ArrayRM2(x));

% h = figure('position',[100, 100, 1000, 750]); set(gcf,'color','w'); set(gca,'FontSize',24);
axes(handles.axes1)

plot(ArrayRT1,ArrayRT2,'ob',ArrayRM1,ArrayRM2,'*r','LineWidth',3),grid on
% axis([0 5000 0 5000])
h_legend = legend('Target trajectory','Missile trajectory','Location','NorthWest');
set(h_legend,'FontSize',14); set(gca,'FontSize',14);
title('Two-dimensional tactical missile-target engagement simulation')
xlabel('Downrange (Ft)','FontSize', 16)
ylabel('Altitude or crossrange (Ft)','FontSize', 16)


% saveTightFigure(h,'trajectory20NN5.pdf')
%===========================================================
% h = figure('position',[100, 100, 1000, 750]); set(gcf,'color','w'); set(gca,'FontSize',24);
% plot(ArrayT,ArrayXNCG,'LineWidth',3),grid on
% %title('Two-dimensional tactical missile-target engagement simulation')
% xlabel('Time (sec)','FontSize', 16)
% ylabel('Acceleration of missle (G)','FontSize', 16)
% saveTightFigure(h,'MissileAcceleration20NN5.pdf')
% output=[ArrayT',ArrayRT1',ArrayRT2',ArrayRM1',ArrayRM2',ArrayXNCG',ArrayRTM' ];
% %save datfil.txt output /ascii
% disp '*** Simulation Complete'
% toc;
hold on
% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%TWO-DIMENSIONAL TACTICAL MISSILE-TARGET ENGAGEMENT SIMULATION
%solving using the second-order Runge–Kutta numerical integration technique

%% Simulation inputs
dx=str2num(get(handles.dx,'string'));
dy=str2num(get(handles.dy,'string'));

RT1=str2num(get(handles.Rt1,'string'))+dx;
RT2=str2num(get(handles.Rt2,'string'));
RM1=str2num(get(handles.Rm1,'string'));
RM2=str2num(get(handles.Rm2,'string'));
VM=str2num(get(handles.VM,'string'));
VT=str2num(get(handles.VT,'string'));
XNP=str2num(get(handles.XNP,'string'));
XNT=str2num(get(handles.XNTT,'string'));
HEDEG=str2num(get(handles.HEDEG,'string'));

set(handles.Rt1,'string',RT1);

n=0; %counter on points
%VM = 3000.; %magnitude of the missile velocity
%VT = 1000.; %magintude of the target velocity
%XNTT = 96.6; %target accelration (target manuver)[first error source]
%HEDEG = 0.0; %heading error (in degrees)       [second error source]
%XNP = 5.; %effictive navigation ratio
%RM1 = 0.; %initial location of the missile (in the dowenrange axis)
%RM2 = 9000.; %initial location of the missile (in the crossrange axis)
%RT1 = 12000.; %initial location of the target (in the dowenrange axis)
%RT2 = 10000.; %initial location of the target (in the crossrange axis)
BETA=0.; %angle between the target velocity vector and the doewrange axis

%% Diffrential equations

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
VC=-(RTM1*VTM1 + RTM2*VTM2)/RTM; %closing velocity

while VC >= 0 %terminate the programme when the velocity chnges its sign
    %means that the separation between the missile and target is a minimum
    if RTM < 1000
        H=.0002; %step size made smaller near the end of the flight
        %(to accurately capture the magnitude of the miss distance)
    else
        H=.01; %step of the most of the flight (except the end)
    end
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
        XNC=XNP*VC*XLAMD;
        XNC1=(XNP+20000000*XLAMD^2/VC)*VC*XLAMD;
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
        ArrayRT1(n)=RT1;
        ArrayRT2(n)=RT2;
        ArrayRM1(n)=RM1;
        ArrayRM2(n)=RM2;
        ArrayXNCG(n)=XNC/32.2;
        ArrayRTM(n)=RTM;
    end
end
RTM %miss distance
x=floor(VM/VT);
axes(handles.axes1)
set(handles.Rm1,'string',ArrayRM1(x));
set(handles.Rm2,'string',ArrayRM2(x));
% h = figure('position',[100, 100, 1000, 750]); set(gcf,'color','w'); set(gca,'FontSize',24);
plot(ArrayRT1,ArrayRT2,'ob',ArrayRM1,ArrayRM2,'*r','LineWidth',3),grid on
% axis([0 5000 0 5000])
h_legend = legend('Target trajectory','Missile trajectory','Location','NorthWest');
set(h_legend,'FontSize',14); set(gca,'FontSize',14);
title('Two-dimensional tactical missile-target engagement simulation')
xlabel('Downrange (Ft)','FontSize', 16)
ylabel('Altitude or crossrange (Ft)','FontSize', 16)
% saveTightFigure(h,'trajectory20NN5.pdf')
%===========================================================
% h = figure('position',[100, 100, 1000, 750]); set(gcf,'color','w'); set(gca,'FontSize',24);
% plot(ArrayT,ArrayXNCG,'LineWidth',3),grid on
%title('Two-dimensional tactical missile-target engagement simulation')
% xlabel('Time (sec)','FontSize', 16)
% ylabel('Acceleration of missle (G)','FontSize', 16)
% saveTightFigure(h,'MissileAcceleration20NN5.pdf')
% output=[ArrayT',ArrayRT1',ArrayRT2',ArrayRM1',ArrayRM2',ArrayXNCG',ArrayRTM' ];
%save datfil.txt output /ascii
% disp '*** Simulation Complete'

% toc;
hold on
% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%TWO-DIMENSIONAL TACTICAL MISSILE-TARGET ENGAGEMENT SIMULATION
%solving using the second-order Runge–Kutta numerical integration technique

%% Simulation inputs
dx=str2num(get(handles.dx,'string'));
dy=str2num(get(handles.dy,'string'));

RT1=str2num(get(handles.Rt1,'string'))-dx;
RT2=str2num(get(handles.Rt2,'string'));
RM1=str2num(get(handles.Rm1,'string'));
RM2=str2num(get(handles.Rm2,'string'));
VM=str2num(get(handles.VM,'string'));
VT=str2num(get(handles.VT,'string'));
XNP=str2num(get(handles.XNP,'string'));
XNT=str2num(get(handles.XNTT,'string'));
HEDEG=str2num(get(handles.HEDEG,'string'));

set(handles.Rt1,'string',RT1);

n=0; %counter on points
%VM = 3000.; %magnitude of the missile velocity
%VT = 1000.; %magintude of the target velocity
%XNTT = 96.6; %target accelration (target manuver)[first error source]
%HEDEG = 0.0; %heading error (in degrees)       [second error source]
%XNP = 5.; %effictive navigation ratio
%RM1 = 0.; %initial location of the missile (in the dowenrange axis)
%RM2 = 9000.; %initial location of the missile (in the crossrange axis)
%RT1 = 12000.; %initial location of the target (in the dowenrange axis)
%RT2 = 10000.; %initial location of the target (in the crossrange axis)
BETA=0.; %angle between the target velocity vector and the doewrange axis

%% Diffrential equations

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
VC=-(RTM1*VTM1 + RTM2*VTM2)/RTM; %closing velocity

while VC >= 0 %terminate the programme when the velocity chnges its sign
    %means that the separation between the missile and target is a minimum
    if RTM < 1000
        H=.0002; %step size made smaller near the end of the flight
        %(to accurately capture the magnitude of the miss distance)
    else
        H=.01; %step of the most of the flight (except the end)
    end
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
        XNC=XNP*VC*XLAMD;
        XNC1=(XNP+20000000*XLAMD^2/VC)*VC*XLAMD;
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
        ArrayRT1(n)=RT1;
        ArrayRT2(n)=RT2;
        ArrayRM1(n)=RM1;
        ArrayRM2(n)=RM2;
        ArrayXNCG(n)=XNC/32.2;
        ArrayRTM(n)=RTM;
    end
end
RTM %miss distance
x=floor(VM/VT);

axes(handles.axes1)
set(handles.Rm1,'string',ArrayRM1);
set(handles.Rm2,'string',ArrayRM2);
% h = figure('position',[100, 100, 1000, 750]); set(gcf,'color','w'); set(gca,'FontSize',24);
plot(ArrayRT1,ArrayRT2,'ob',ArrayRM1,ArrayRM2,'*r','LineWidth',3),grid on
% axis([0 5000 0 5000])
h_legend = legend('Target trajectory','Missile trajectory','Location','NorthWest');
set(h_legend,'FontSize',14); set(gca,'FontSize',14);
title('Two-dimensional tactical missile-target engagement simulation')
xlabel('Downrange (Ft)','FontSize', 16)
ylabel('Altitude or crossrange (Ft)','FontSize', 16)
% saveTightFigure(h,'trajectory20NN5.pdf')
%===========================================================
% h = figure('position',[100, 100, 1000, 750]); set(gcf,'color','w'); set(gca,'FontSize',24);
% plot(ArrayT,ArrayXNCG,'LineWidth',3),grid on
%title('Two-dimensional tactical missile-target engagement simulation')
% xlabel('Time (sec)','FontSize', 16)
% ylabel('Acceleration of missle (G)','FontSize', 16)
% saveTightFigure(h,'MissileAcceleration20NN5.pdf')
% output=[ArrayT',ArrayRT1',ArrayRT2',ArrayRM1',ArrayRM2',ArrayXNCG',ArrayRTM' ];
%save datfil.txt output /ascii
% disp '*** Simulation Complete'

% toc;

% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%TWO-DIMENSIONAL TACTICAL MISSILE-TARGET ENGAGEMENT SIMULATION
%solving using the second-order Runge–Kutta numerical integration technique

%% Simulation inputs
dx=str2num(get(handles.dx,'string'));
dy=str2num(get(handles.dy,'string'));

RT1=str2num(get(handles.Rt1,'string'));
RT2=str2num(get(handles.Rt2,'string'))-dy;
RM1=str2num(get(handles.Rm1,'string'));
RM2=str2num(get(handles.Rm2,'string'));
VM=str2num(get(handles.VM,'string'));
VT=str2num(get(handles.VT,'string'));
XNP=str2num(get(handles.XNP,'string'));
XNT=str2num(get(handles.XNTT,'string'));
HEDEG=str2num(get(handles.HEDEG,'string'));

set(handles.Rt2,'string',RT2);

n=0; %counter on points
%VM = 3000.; %magnitude of the missile velocity
%VT = 1000.; %magintude of the target velocity
%XNTT = 96.6; %target accelration (target manuver)[first error source]
%HEDEG = 0.0; %heading error (in degrees)       [second error source]
%XNP = 5.; %effictive navigation ratio
%RM1 = 0.; %initial location of the missile (in the dowenrange axis)
%RM2 = 9000.; %initial location of the missile (in the crossrange axis)
%RT1 = 12000.; %initial location of the target (in the dowenrange axis)
%RT2 = 10000.; %initial location of the target (in the crossrange axis)
BETA=0.; %angle between the target velocity vector and the doewrange axis

%% Diffrential equations

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
VC=-(RTM1*VTM1 + RTM2*VTM2)/RTM; %closing velocity

while VC >= 0 %terminate the programme when the velocity chnges its sign
    %means that the separation between the missile and target is a minimum
    if RTM < 1000
        H=.0002; %step size made smaller near the end of the flight
        %(to accurately capture the magnitude of the miss distance)
    else
        H=.01; %step of the most of the flight (except the end)
    end
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
        XNC=XNP*VC*XLAMD;
        XNC1=(XNP+20000000*XLAMD^2/VC)*VC*XLAMD;
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
        ArrayRT1(n)=RT1;
        ArrayRT2(n)=RT2;
        ArrayRM1(n)=RM1;
        ArrayRM2(n)=RM2;
        ArrayXNCG(n)=XNC/32.2;
        ArrayRTM(n)=RTM;
    end
end
RTM %miss distance
x=floor(VM/VT);
axes(handles.axes1)
set(handles.Rm1,'string',ArrayRM1(x));
set(handles.Rm2,'string',ArrayRM2(x));
% h = figure('position',[100, 100, 1000, 750]); set(gcf,'color','w'); set(gca,'FontSize',24);
plot(ArrayRT1,ArrayRT2,'ob',ArrayRM1,ArrayRM2,'*r','LineWidth',3),grid on
% axis([0 5000 0 5000])
h_legend = legend('Target trajectory','Missile trajectory','Location','NorthWest');
set(h_legend,'FontSize',14); set(gca,'FontSize',14);
title('Two-dimensional tactical missile-target engagement simulation')
xlabel('Downrange (Ft)','FontSize', 16)
ylabel('Altitude or crossrange (Ft)','FontSize', 16)
% saveTightFigure(h,'trajectory20NN5.pdf')
%===========================================================
% h = figure('position',[100, 100, 1000, 750]); set(gcf,'color','w'); set(gca,'FontSize',24);
% plot(ArrayT,ArrayXNCG,'LineWidth',3),grid on
%title('Two-dimensional tactical missile-target engagement simulation')
% xlabel('Time (sec)','FontSize', 16)
% ylabel('Acceleration of missle (G)','FontSize', 16)
% saveTightFigure(h,'MissileAcceleration20NN5.pdf')
% output=[ArrayT',ArrayRT1',ArrayRT2',ArrayRM1',ArrayRM2',ArrayXNCG',ArrayRTM' ];
%save datfil.txt output /ascii
% disp '*** Simulation Complete'

% toc;


function XNTT_Callback(hObject, eventdata, handles)
% hObject    handle to XNTT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of XNTT as text
%        str2double(get(hObject,'String')) returns contents of XNTT as a double


% --- Executes during object creation, after setting all properties.
function XNTT_CreateFcn(hObject, eventdata, handles)
% hObject    handle to XNTT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dx_Callback(hObject, eventdata, handles)
% hObject    handle to dx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dx as text
%        str2double(get(hObject,'String')) returns contents of dx as a double


% --- Executes during object creation, after setting all properties.
function dx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dy_Callback(hObject, eventdata, handles)
% hObject    handle to dy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dy as text
%        str2double(get(hObject,'String')) returns contents of dy as a double


% --- Executes during object creation, after setting all properties.
function dy_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function n_Callback(hObject, eventdata, handles)
% hObject    handle to n (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of n as text
%        str2double(get(hObject,'String')) returns contents of n as a double


% --- Executes during object creation, after setting all properties.
function n_CreateFcn(hObject, eventdata, handles)
% hObject    handle to n (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.Rt1,'string',40000)
set(handles.Rt2,'string',10000)
set(handles.Rm1,'string',1)
set(handles.Rm2,'string',1)
set(handles.VM,'string',3000)
set(handles.VT,'string',1000)
set(handles.HEDEG,'string',0)
set(handles.XNP,'string',3)
set(handles.XNTT,'string',5)

set(handles.N,'string',40)
set(handles.W1,'string',1)
set(handles.W2,'string',100000)
set(handles.degree,'string',4)
set(handles.amax1,'string',12)
set(handles.t0,'string',0)
set(handles.t3,'string',10)


axes(handles.axes1)
plot(0,0)
axes(handles.axes2)
plot(0,0)
axes(handles.axes3)
plot(0,0)


function times_Callback(hObject, eventdata, handles)
% hObject    handle to times (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of times as text
%        str2double(get(hObject,'String')) returns contents of times as a double


% --- Executes during object creation, after setting all properties.
function times_CreateFcn(hObject, eventdata, handles)
% hObject    handle to times (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function pushbutton2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on selection change in pop1.
function pop1_Callback(hObject, eventdata, handles)
% hObject    handle to pop1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pop1=get(handles.pop1,'value');
w1=str2num(get(handles.W1,'string'));
w2=str2num(get(handles.W2,'string'));
N=str2num(get(handles.N,'string'));

if pop1==2
    set(handles.pop2,'visible','on');
    set(handles.text22,'visible','on');
    set(handles.N,'visible','on');
    set(handles.W1,'visible','on');
    set(handles.W2,'visible','on');
    set(handles.text23,'visible','on');
    set(handles.text24,'visible','on');
    set(handles.text25,'visible','on');
elseif pop1==1
    set(handles.pop2,'visible','off');
    set(handles.text22,'visible','off');
    set(handles.N,'visible','off');
    set(handles.W1,'visible','off');
    set(handles.W2,'visible','off');
    set(handles.text23,'visible','off');
    set(handles.text24,'visible','off');
    set(handles.text25,'visible','off');
    set(handles.text26,'visible','off');
    set(handles.text27,'visible','off');
    set(handles.text28,'visible','off');
    set(handles.text29,'visible','off');
    set(handles.degree,'visible','off');
    set(handles.amax1,'visible','off');
    set(handles.t0,'visible','off');
    set(handles.t3,'visible','off');
end
% Hints: contents = cellstr(get(hObject,'String')) returns pop1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop1

% --- Executes during object creation, after setting all properties.
function pop1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on selection change in pop2.
function pop2_Callback(hObject, eventdata, handles)
% hObject    handle to pop2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pop2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop2
pop2=get(handles.pop2,'value');
if pop2==2
    set(handles.degree,'visible','on');
    set(handles.text26,'visible','on');
    
elseif pop2==1
    set(handles.amax1,'visible','off');
    set(handles.text27,'visible','off');
    set(handles.t0,'visible','off');
    set(handles.text28,'visible','off');
    set(handles.t3,'visible','off');
    set(handles.text29,'visible','off');
    set(handles.degree,'visible','off');
    set(handles.text26,'visible','off');
elseif pop2==3
    set(handles.amax1,'visible','on');
    set(handles.text27,'visible','on');
    set(handles.t0,'visible','on');
    set(handles.text28,'visible','on');
    set(handles.t3,'visible','on');
    set(handles.text29,'visible','on');
    set(handles.degree,'visible','off');
    set(handles.text26,'visible','off');
elseif pop2==4
    set(handles.amax1,'visible','on');
    set(handles.text27,'visible','on');
    set(handles.t0,'visible','on');
    set(handles.text28,'visible','on');
    set(handles.t3,'visible','on');
    set(handles.text29,'visible','on');
    set(handles.degree,'visible','off');
    set(handles.text26,'visible','off');
end
% --- Executes during object creation, after setting all properties.
function pop2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on key press with focus on pop1 and none of its controls.
function pop1_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to pop1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
pop1=get(handles.pop1,'value');

if pop1==2
    set(handles.pop2,'visible','on');
end


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over pop1.
function pop1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to pop1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pop1=get(handles.pop1,'value');

if pop1==2
    set(handles.pop2,'visible','on');
end


% --- Executes during object creation, after setting all properties.
function run_CreateFcn(hObject, eventdata, handles)
% hObject    handle to run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

function [a] =  trapacc(~,t1,t2,t3,amax,T)

%%Parameters of the trapezoid function
%amax1=max accelartion
%t0=beginning of the ramp
%t1=end of the ramp & beginning of the const. acc.
%t2=end of the const. acc. & beginning of the decent
%t3=end of the decent

T=mod(T,2*t3);

if T<t1
    a= (amax/t1)*T;
elseif T<t2
    a= amax;
elseif T<2*t3-t2
    a= (amax/(t2-t3))*(T-t3);
elseif T<2*t3-t1
    a= -amax;
elseif T<2*t3
    a= amax *(T-2*t3);
    
elseif 2*t3<T
    error('time value out of range');
    
end




function [RTM,XNCmax,ArrayRT1b,ArrayRT2b,ArrayRM1b,ArrayRM2b,ArrayTb,ArrayXNCGb,ArrayXNTGb] = PNeq (polyORtrapORtrapSymm,polydegree, c , draw);

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
RT1 = 40000.; %initial location of the target (in the dowenrange axis)[ft]
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
VC=-(RTM1*VTM1 + RTM2*VTM2)/RTM; %closing velocity

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
ArrayRT1b=ArrayRT1;
ArrayRT2b=ArrayRT2;
ArrayRM1b=ArrayRM1;
ArrayRM2b=ArrayRM2;
ArrayTb=ArrayT;
ArrayXNCGb=ArrayXNCG;
ArrayXNTGb=ArrayXNTG;

%%----outputs----
% if draw==1 %to viwe all figurs (runs)
%     h = figure('position',[100, 100, 1000, 750]); set(gcf,'color','w'); set(gca,'FontSize',24);
%
%
% plot(ArrayRT1,ArrayRT2,':',ArrayRM1,ArrayRM2,'*r','LineWidth',3),grid on
%     h_legend = legend('Target trajectory','Missile trajectory','Location','NorthWest');
%     set(h_legend,'FontSize',14); set(gca,'FontSize',14);
%     title('Two-dimensional tactical missile-target engagement simulation')
%     xlabel('Downrange (Ft) ','FontSize', 16)
%     ylabel('Altitude or crossrange (Ft)','FontSize', 16)
% %     saveTightFigure(h,'trajectoryTRAP2.pdf')
%
%     hold on   %-------------------------------
%
%     h = figure('position',[100, 100, 1000, 750]); set(gcf,'color','w'); set(gca,'FontSize',24);
% % axes(handles.axes2)
% plot(ArrayT,ArrayXNCG,'LineWidth',3),grid on
%     %title('Two-dimensional tactical missile-target engagement simulation')
%     xlabel('Time (sec)','FontSize', 16)
%     ylabel('Acceleration of missle (G)','FontSize', 16)
% %    saveTightFigure(h,'MissileAccelerationTRAP2.pdf')
%
%     hold on
%
%     h = figure('position',[100, 100, 1000, 750]); set(gcf,'color','w'); set(gca,'FontSize',24);
% % axes(handles.axes3)
% plot(ArrayT,ArrayXNTG,'LineWidth',3),grid on
%     %title('Two-dimensional tactical missile-target engagement simulation')
%     xlabel('Time (sec)','FontSize', 16)
%     ylabel('Acceleration of Target (G)','FontSize', 16)
% %     saveTightFigure(h,'TargetAccelerationTRAP2.pdf')
% end

% --- Executes on button press in run.
function run_Callback(hObject, eventdata, handles)
% hObject    handle to run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% set(handles.push,'string',push)
% set(handles.time,'string',t)
hold off

G=32.2;
pop1=get(handles.pop1,'value');
pop2=get(handles.pop2,'value');
w1=str2num(get(handles.W1,'string'));
w2=str2num(get(handles.W2,'string'));
N=str2num(get(handles.N,'string'));
amax11=str2num(get(handles.amax1,'string'));
amax=amax11*G;
t0=str2num(get(handles.t0,'string'));
t3=str2num(get(handles.t3,'string'));

if pop1==2   %drop down menu for pn optimal
    
    dx=str2num(get(handles.dx,'string'));
    dy=str2num(get(handles.dy,'string'));
    
    RT1=str2num(get(handles.Rt1,'string'));
    RT2=str2num(get(handles.Rt2,'string'));
    RM1=str2num(get(handles.Rm1,'string'));
    RM2=str2num(get(handles.Rm2,'string'));
    VM=str2num(get(handles.VM,'string'));
    VT=str2num(get(handles.VT,'string'));
    XNP=str2num(get(handles.XNP,'string'));
    XNTT=str2num(get(handles.XNTT,'string'));
    HEDEG=str2num(get(handles.HEDEG,'string'));
    
    set(handles.Rt2,'string',RT2);
    
    % set(handles.axes2,'visible','on');
    
    
    n=0; %counter on points
    G=32.2;
    % VM = 3000.; %magnitude of the missile velocity
    % VT = 1000.; %magintude of the target velocity
    XNT =XNTT*G ; %target accelration (target manuver)[first error source]
    % HEDEG = 0.0; %heading error (in degrees)       [second error source]
    % XNP = 5.; %effictive navigation ratio
    % RM1 = 0.; %initial location of the missile (in the dowenrange axis)
    % RM2 = 0.; %initial location of the missile (in the crossrange axis)
    % RT1 = 40000.; %initial location of the target (in the dowenrange axis)
    % RT2 = 10000.; %initial location of the target (in the crossrange axis)
    BETA=0.; %angle between the target velocity vector and the doewrange axis
    
    %% Diffrential equations
    
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
    VC=-(RTM1*VTM1 + RTM2*VTM2)/RTM; %closing velocity
    
    while VC >= 0 %terminate the programme when the velocity chnges its sign
        %means that the separation between the missile and target is a minimum
        if RTM < 1000
            H=.0002; %step size made smaller near the end of the flight
            %(to accurately capture the magnitude of the miss distance)
        else
            H=.01; %step of the most of the flight (except the end)
        end
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
            XNC=XNP*VC*XLAMD;
            %XNC1=(XNP+20000000*XLAMD^2/VC)*VC*XLAMD;
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
            ArrayXNCG(n)=XNC/G;
            ArrayRTM(n)=RTM;
        end
    end
    set(handles.RTMM,'string',RTM) %miss distance
    axes(handles.axes1)
    plot(ArrayRT1,ArrayRT2,':',ArrayRM1,ArrayRM2,'*r','LineWidth',3),grid on
    % figure('position',[100, 100, 1000, 750]); set(gcf,'color','w'); set(gca,'FontSize',24);
    legend('Target trajectory','Missile trajectory','Location','NorthWest');
    % set(h_legend,'FontSize',14); set(gca,'FontSize',14);
    title('Two-dimensional tactical missile-target engagement simulation')
    xlabel('Downrange (Ft)','FontSize', 16)
    ylabel('Altitude or crossrange (Ft)','FontSize', 16)
    % saveTightFigure(h,'trajectoryXNT5HE0N5.pdf')
    %===========================================================
    axes(handles.axes2)
    % h = figure('position',[100, 100, 1000, 750]); set(gcf,'color','w'); set(gca,'FontSize',24);
    plot(ArrayT,ArrayXNCG,'LineWidth',3),grid on
    %title('Two-dimensional tactical missile-target engagement simulation')
    xlabel('Time (sec)','FontSize', 16)
    ylabel('Acceleration of missle (G)','FontSize', 16)
    % saveTightFigure(h,'MissileAccelerationXNT5HE0N5.pdf')
    output=[ArrayT',ArrayRT1',ArrayRT2',ArrayRM1',ArrayRM2',ArrayXNCG',ArrayRTM' ];
    %save datfil.txt output /ascii
    %==========================================================
    axes(handles.axes3)
    % h = figure('position',[100, 100, 1000, 750]); set(gcf,'color','w'); set(gca,'FontSize',24);
    plot(ArrayT,ArrayXNTG,'LineWidth',3),grid on
    %title('Two-dimensional tactical missile-target engagement simulation')
    xlabel('Time (sec)','FontSize', 16)
    ylabel('Acceleration of Target (G)','FontSize', 16)
    % saveTightFigure(h,'TargetAccelerationXNT5HE0N5.pdf')
    
    % disp '*** Simulation Complete'
    
    
    
    %%%%%%%%%##########$$$$$$$$$$$$$$$$$&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    
    %This code used PN eq to get 3 types of evasive maneuvers :
    % First case [poly] :  we make the escape maneuver as nth degree polynomial
    % and we try to find the optimized solution using PN eq
    % Second case [Trap-general] :  we make the escape maneuver as trapezoidal
    % acceleration, and we sacn all the domain for t1 and t2 which lead to
    % otimized solution
    % Third case [Tooth - square] : special case of trap - symmetric
    %% Configuration Parameters
    
    pop2=get(handles.pop2,'value');
    
    if pop2==2
        polyORtrapORtrapSymm = 'poly'; % choose which type of target acceleration
    elseif pop2==3
        polyORtrapORtrapSymm = 'trap'; % choose which type of target acceleration
    elseif pop2==4
        polyORtrapORtrapSymm = 'trapSymm'; % choose which type of target acceleration
    end
    
    
    % N=20; %number of runs
    G = 32.2; %gravity
    
    
    %=====Parameters of the polynomial function=====
    degree = 5; %Polynomial degree
    
    %=====Parameters of the trapezoid function=====
    % amax1=12*G; %max accelartion
    % t0=0;  %beginning of the ramp
    %t1  %end of the ramp & beginning of the const. acc.
    %t2  %end of the const. acc. & beginning of the decent
    % t3=6;  %end of the decent
    
    %---cost function weighting---
    % f=w1*XNC + w2*RTM;
    
    
    %% Prepareing the variables for each case :
    
    if strcmp(polyORtrapORtrapSymm,'poly')
        cs=zeros(N,degree); % Matrix of Coefficients: generating vector of zeros
        f=zeros(N,1); % cost array: generating vector of zeros
        
    elseif strcmp(polyORtrapORtrapSymm,'trap')
        cs=zeros(N,N,5); % Matrix of coefficients trap ; 5 is length(c)
        f=zeros(N,N); % cost matrix: generating matrix of zeros
        
        t1 = zeros(N,1);
        t2 = zeros(N,N);
    elseif strcmp(polyORtrapORtrapSymm,'trapSymm')
        cs=zeros(N,5); % Matrix of coefficients trap ; 5 is length(c)
        f=zeros(N,1); % cost array: generating vector of zeros
        
        t1 = zeros(N,1);
        t2 = zeros(N,1);
    end
    
    %% Three cases of evasive maneuvers
    
    if strcmp(polyORtrapORtrapSymm,'poly') %===== [1] POLY SOLUTION =====
        
        for(i=1:1:N)
            %making the escape manuver for the target as a nth degree polynomial
            c = 3* G* rand(1,degree); %random generation for the coeff. of the polynomial
            
            cs(i,:)=c;
            [RTM,XNCmax,ArrayRT1b,ArrayRT2b,ArrayRM1b,ArrayRM2b,ArrayTb,ArrayXNCGb,ArrayXNTGb] = PNeq (polyORtrapORtrapSymm, degree, c ,0);
            
            %---cost function---
            % filling the vector with the values of the cost fn
            f(i,1)=w1*(XNCmax^2) + w2*(RTM^2)
            
        end
        
        [s,R]= max(f) %get index and value of max. f
        cs(R,:); %the numbers that generate max cost
        [RTM,XNCmax,ArrayRT1b,ArrayRT2b,ArrayRM1b,ArrayRM2b,ArrayTb,ArrayXNCGb,ArrayXNTGb] = PNeq (polyORtrapORtrapSymm, degree, cs(R,:),1); %draw the optimized solution
        
        %-------------------------------------------------------------------------
        
    elseif strcmp(polyORtrapORtrapSymm,'trap')  %=== [2] TRAPEZOIDAL SOLUTION ===
        %         c= rand(1,2); %generating array of two rundom number
        %         ct=sort(c); %rearange the array, starting from the smallest value
        %         t1= t3*ct(1); %t1>t2 & scaling with t3
        %         t2 = t3*ct(2);
        
        for(i=1:1:N+1)
            t1(i)= t0 + i*((t3-t0)/N);
            
            for (j=i:1:N+1)
                t2(i,j)= t1(i,1)+ j*((t3-t0)/N);
                c = [t0 t1(i,1) t2(i,j) t3 amax];
                cs(i,j,:)= c;
                [RTM,XNCmax,ArrayRT1b,ArrayRT2b,ArrayRM1b,ArrayRM2b,ArrayTb,ArrayXNCGb,ArrayXNTGb] = PNeq (polyORtrapORtrapSymm, degree, c ,0);
                
                %---cost function---
                % filling the vector with the values of the cost fn
                f(i,j)= w1*(XNCmax^2) + w2*(RTM^2);
            end
            
        end
        
        %f(:); %put all the matrix elements in column array order
        [M,I] = max(f(:)); %get the max of this array and its index
        [I_row, I_col] = ind2sub(size(f),I); %refer the index to the original matrix
        
        %[s,R]= max(f)
        
        cs(I_row, I_col,:); %the numbers that generate max cost
        [RTM,XNCmax,ArrayRT1b,ArrayRT2b,ArrayRM1b,ArrayRM2b,ArrayTb,ArrayXNCGb,ArrayXNTGb] = PNeq (polyORtrapORtrapSymm, degree, cs(I_row, I_col,:),1);
        %I_row
        %I_col
        
        %-------------------------------------------------------------------
        
    elseif strcmp(polyORtrapORtrapSymm,'trapSymm')  %=== [3] Tooth - Square SOLUTION (speacial case of trap)===
        
        for(i=1:1:N+1/2)
            t1(i)= t0 + (i-1)*((t3-t0)/N);
            t2(i)= t3 - (i-1)*((t3-t0)/N);
            
            c = [t0 t1(i) t2(i) t3 amax];
            cs(i,:)= c;
            [RTM,XNCmax,ArrayRT1b,ArrayRT2b,ArrayRM1b,ArrayRM2b,ArrayTb,ArrayXNCGb,ArrayXNTGb] = PNeq (polyORtrapORtrapSymm, degree, c ,0);
            
            %---cost function---
            % filling the vector with the values of the cost fn
            f(i,1)= w1*(XNCmax^2) + w2*(RTM^2);
            
            
        end
        
        [s,R]= max(f) %get index and value of max. f
        cs(R,:); %the numbers that generate max cost
        [RTM,XNCmax,ArrayRT1b,ArrayRT2b,ArrayRM1b,ArrayRM2b,ArrayTb,ArrayXNCGb,ArrayXNTGb] = PNeq (polyORtrapORtrapSymm, degree, cs(R,:),1) %draw the optimized solution
    end
    
    %output=[ArrayT',ArrayRT1',ArrayRT2',ArrayRM1',ArrayRM2',ArrayXNCG',ArrayRTM' ];
    %save datfil.txt output /ascii
    % disp '*** Simulation Complete'
    axes(handles.axes1)
    plot(ArrayRT1b,ArrayRT2b,':',ArrayRM1b,ArrayRM2b,'*r','LineWidth',3),grid on
    h_legend = legend('Target trajectory','Missile trajectory','Location','NorthWest');
    set(h_legend,'FontSize',10); set(gca,'FontSize',10);
    title('Two-dimensional tactical missile-target engagement simulation')
    xlabel('Downrange (Ft) ','FontSize', 10)
    ylabel('Altitude or crossrange (Ft)','FontSize', 10)
    
    axes(handles.axes2)
    plot(ArrayTb,ArrayXNCGb,'LineWidth',3),grid on
    %title('Two-dimensional tactical missile-target engagement simulation')
    xlabel('Time (sec)','FontSize', 10)
    ylabel('Acceleration of missle (G)','FontSize', 10)
    
    axes(handles.axes3)
    plot(ArrayTb,ArrayXNTGb,'LineWidth',3),grid on
    %     title('Two-dimensional tactical missile-target engagement simulation')
    xlabel('Time (sec)','FontSize', 10)
    ylabel('Acceleration of Target (G)','FontSize', 10)
    
elseif pop1==1
    set(handles.axes2,'visible','off')
    set(handles.axes3,'visible','off')
elseif pop1==3     % optimal
    set(handles.axes2,'visible','off')
    set(handles.axes3,'visible','off')
    
elseif pop1==4    % Augmented
    set(handles.axes2,'visible','off')
    set(handles.axes3,'visible','off')
end



function N_Callback(hObject, eventdata, handles)
% hObject    handle to N (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of N as text
%        str2double(get(hObject,'String')) returns contents of N as a double


% --- Executes during object creation, after setting all properties.
function N_CreateFcn(hObject, eventdata, handles)
% hObject    handle to N (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function W1_Callback(hObject, eventdata, handles)
% hObject    handle to W1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of W1 as text
%        str2double(get(hObject,'String')) returns contents of W1 as a double


% --- Executes during object creation, after setting all properties.
function W1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to W1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function W2_Callback(hObject, eventdata, handles)
% hObject    handle to W2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of W2 as text
%        str2double(get(hObject,'String')) returns contents of W2 as a double


% --- Executes during object creation, after setting all properties.
function W2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to W2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function degree_Callback(hObject, eventdata, handles)
% hObject    handle to degree (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of degree as text
%        str2double(get(hObject,'String')) returns contents of degree as a double


% --- Executes during object creation, after setting all properties.
function degree_CreateFcn(hObject, eventdata, handles)
% hObject    handle to degree (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function amax1_Callback(hObject, eventdata, handles)
% hObject    handle to amax1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of amax1 as text
%        str2double(get(hObject,'String')) returns contents of amax1 as a double


% --- Executes during object creation, after setting all properties.
function amax1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to amax1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function t0_Callback(hObject, eventdata, handles)
% hObject    handle to t0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of t0 as text
%        str2double(get(hObject,'String')) returns contents of t0 as a double


% --- Executes during object creation, after setting all properties.
function t0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to t0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function t3_Callback(hObject, eventdata, handles)
% hObject    handle to t3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of t3 as text
%        str2double(get(hObject,'String')) returns contents of t3 as a double


% --- Executes during object creation, after setting all properties.
function t3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to t3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
