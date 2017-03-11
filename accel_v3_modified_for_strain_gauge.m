function varargout = accel_v3_modified_for_strain_gauge(varargin)
% ACCEL_V3_MODIFIED_FOR_STRAIN_GAUGE MATLAB code for accel_v3_modified_for_strain_gauge.fig
%      ACCEL_V3_MODIFIED_FOR_STRAIN_GAUGE, by itself, creates a new ACCEL_V3_MODIFIED_FOR_STRAIN_GAUGE or raises the existing
%      singleton*.
%
%      H = ACCEL_V3_MODIFIED_FOR_STRAIN_GAUGE returns the handle to a new ACCEL_V3_MODIFIED_FOR_STRAIN_GAUGE or the handle to
%      the existing singleton*.
%
%      ACCEL_V3_MODIFIED_FOR_STRAIN_GAUGE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ACCEL_V3_MODIFIED_FOR_STRAIN_GAUGE.M with the given input arguments.
%
%      ACCEL_V3_MODIFIED_FOR_STRAIN_GAUGE('Property','Value',...) creates a new ACCEL_V3_MODIFIED_FOR_STRAIN_GAUGE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before accel_v3_modified_for_strain_gauge_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to accel_v3_modified_for_strain_gauge_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help accel_v3_modified_for_strain_gauge

% Last Modified by GUIDE v2.5 18-Nov-2016 16:04:30

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @accel_v3_modified_for_strain_gauge_OpeningFcn, ...
    'gui_OutputFcn',  @accel_v3_modified_for_strain_gauge_OutputFcn, ...
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


% --- Executes just before accel_v3_modified_for_strain_gauge is made visible.
function accel_v3_modified_for_strain_gauge_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to accel_v3_modified_for_strain_gauge (see VARARGIN)

% Choose default command line output for accel_v3_modified_for_strain_gauge
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes accel_v3_modified_for_strain_gauge wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = accel_v3_modified_for_strain_gauge_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function clearUserData(handles)
set(handles.Import_accel,'UserData',[]);
set(handles.Import_c40,'UserData',[],'Enable','on');
set(handles.Freq_textbox,'UserData',[]);
set(handles.Align_encoder,'UserData',[],'Enable','on');
set(handles.ReAlign,'UserData',[],'Enable','on');
set(handles.Fine_tune_encoder,'UserData',...
    get(handles.Fine_tune_encoder,'Value'),'Enable','on');
set(handles.Plot,'UserData',[]);
set(handles.Calc_vel,'UserData',[]);
set(handles.eventlist,'UserData',[]);
set(handles.AccelCh,'UserData',[]);
set(handles.Zoom,'UserData',[]);
set(handles.axes2,'XScale','linear','YScale','linear');
cla(handles.axes1);
cla(handles.axes2);

% --- Executes on button press in Import_accel.
function Import_accel_Callback(hObject, eventdata, handles)
% hObject    handle to Import_accel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clearUserData(handles);
eventlist_Callback(handles.eventlist, eventdata, handles);

[AccelFile,AccelPath] = uigetfile(...
    'C:\Users\User\Documents\Test files results on fast cards\*.*','Select high freq file');
if AccelFile == 0
    return;
end
set(hObject,'TooltipString',[AccelPath,AccelFile]);
accelFile.AccelFile=AccelFile;
accelFile.AccelPath=AccelPath;
set(hObject,'UserData',accelFile);
guidata(hObject,handles);
Import_c40_Callback(handles.Import_c40, eventdata, handles);


% --- Executes on button press in Import_c40.
function Import_c40_Callback(hObject, eventdata, handles)
% hObject    handle to Import_c40 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[ROGAFileName,ROGAPath] = uigetfile(...
    'C:\Users\User\Documents\Copy of test files on ROGA1\*.c40','Select low freq file to match time');
if ROGAPath == 0
    return;
end
set(hObject,'TooltipString',[ROGAPath,ROGAFileName]);
ROGAFile.ROGAFileName = ROGAFileName;
ROGAFile.ROGAPath = ROGAPath;
set(hObject,'UserData',ROGAFile);

fileID = fopen([ROGAPath,ROGAFileName]);
delimiter = '\t';
formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f';
ROGAc40cell = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines', 1, 'ReturnOnError', false);
fclose(fileID);

ROGATime_s = ROGAc40cell{1,21};
ROGAFreq = 1/ROGATime_s(2,1);
set(handles.Freq_textbox,'UserData',ROGAFreq);  % store low freq in the accel freq textbox just for convenience
ROGAc40matwithNaN = cell2mat(ROGAc40cell);
ROGAc40mat = ROGAc40matwithNaN(3:end-2,:);
set(handles.Plot,'UserData',ROGAc40mat);
guidata(hObject,handles);
Align_encoder_Callback(handles.Align_encoder, eventdata, handles);


function Freq_textbox_Callback(hObject, eventdata, handles)
% hObject    handle to Freq_textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Freq_textbox as text
%        str2double(get(hObject,'String')) returns contents of Freq_textbox as a double


% --- Executes during object creation, after setting all properties.
function Freq_textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Freq_textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Switch_res.
function Switch_res_Callback(hObject, eventdata, handles)
% hObject    handle to Switch_res (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Plot_Callback(handles.Plot, eventdata, handles);
% Hint: get(hObject,'Value') returns toggle state of Switch_res

function Ysn = getYsn(handles,eventdata)
ydatalist_Callback(handles.ydatalist, eventdata, handles);
ydatalist2_Callback(handles.ydatalist2, eventdata, handles);
ydatalist3_Callback(handles.ydatalist3, eventdata, handles);
ydatalist4_Callback(handles.ydatalist4, eventdata, handles);
Ysn = ones(1,4);
Ysn(1) = get(handles.ydatalist,'UserData');
Ysn(2) = get(handles.ydatalist2,'UserData');
Ysn(3) = get(handles.ydatalist3,'UserData');
Ysn(4) = get(handles.ydatalist4,'UserData');

% --- Executes on button press in Plot.
function Plot_Callback(hObject, eventdata, handles)
% hObject    handle to Plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
DB = get(handles.eventlist,'UserData');
if ~isempty(DB)
    set(handles.Fine_tune_encoder,'Enable','off');
    set(handles.ReAlign,'Enable','off');
    set(handles.Align_encoder,'Enable','off');
    set(handles.Import_c40,'Enable','off');
end

plotmethod_Callback(handles.plotmethod, eventdata, handles);

Accelrawdatamat = get(handles.Calc_vel,'UserData');
if isempty(Accelrawdatamat); return; end
AccelOnly = 0;
ROGAc40mat = get(hObject,'UserData');
if isempty(ROGAc40mat); AccelOnly = 1; end
fullres = get(handles.Full_res,'Value');
color = {'b','r','g','m'};
Ysn = getYsn(handles,eventdata);
chop = get(handles.chop,'Value');
AccelFreq = str2double(get(handles.Freq_textbox,'String'));
if isnan(AccelFreq)
    error('Frequency invalid.');
end
try
    offsetSlope = eval(get(handles.offsetAndSlope,'String'));
catch ME
    msgbox('Invalid offset and slope format!');
    return;
end
if ~isnumeric(offsetSlope) || size(offsetSlope,1) ~= 4 || size(offsetSlope,2) ~= 2
    msgbox('Not numeric or size wrong.');
    return;
end
xl = get(handles.axes1,'XLim');
showCh = getShowCh(handles);

cla(handles.axes1);
for i=1:4   % four lines to display
    if Ysn(i) == 0 || showCh(i) == 0
        continue;
    end
    hold(handles.axes1,'on');
    if Ysn(i) > 100
        
        ylistsn = Ysn(i) - 100;
        if ylistsn < 21
            xlistsn = 5*ceil(ylistsn/5)-4;
        else
            xlistsn = 1;
        end
        if fullres == 1
            xarray = Accelrawdatamat(1:end-AccelFreq*chop,xlistsn);
            yarray = (Accelrawdatamat(1:end-AccelFreq*chop,ylistsn)-offsetSlope(i,1)) ./ offsetSlope(i,2);
        else
            xarray = downsample(...
                Accelrawdatamat(1:end-AccelFreq*chop,xlistsn),10);
            yarray = (downsample(...
                Accelrawdatamat(1:end-AccelFreq*chop,ylistsn),10)-offsetSlope(i,1)) ./ offsetSlope(i,2);
        end
        if get(handles.Smooth,'Value') == get(handles.Smooth,'Max')
            yarray = medfilt2(yarray,[3 1]);
        end
    else
        if AccelOnly == 1
            msgbox('ROGA data not available. Please select only Accel channels.');
            return
        end
        ylistsn = Ysn(i) - 0;
        xlistsn = 21;
        xarray = ROGAc40mat(:,xlistsn);
        yarray = (ROGAc40mat(:,ylistsn)-offsetSlope(i,1)) ./ offsetSlope(i,2);
    end
    
    if get(handles.Scaled,'Value') == 1
        ind1 = find(xarray > xl(1),1);
        ind2 = find(xarray > xl(2),1);
        if isempty(ind2) || isempty(ind1); ind1 = 1; ind2 = length(xarray); end
        plot(handles.axes1,xarray,...
            (yarray-mean(yarray(ind1:ind2)))./max(abs(yarray(ind1:ind2)-mean(yarray(ind1:ind2))))...
            ,color{i});
    else
        plot(handles.axes1,xarray,yarray,color{i});
    end
    hold(handles.axes1,'off');
end
set(handles.axes1,'YLimMode','auto');
zoom reset;
guidata(hObject,handles);


function encoderjumpavg = getCardsTimeDiff(test)
% parameters to change before executing
encoderShift = 3;
encoderSpacingThrsh = 5e4;
encoderXcorrWindow = 5e5;

% execution begins here
% calc timeshift of 4 cards
encodercurrent = 1;
encoderjump = 1;
data1 = test{1};
encoderjumpavg = zeros(1,length(test) - 1);
% get to sliding period (encoder-tight zone)
for j = 2:length(data1)
    if (data1(j-1,1) - encoderShift)*...
            (data1(j,1) - encoderShift) < 0
        encoderbefore = encodercurrent;
        encodercurrent = j;
        spacing = encodercurrent - encoderbefore;
        if spacing < encoderSpacingThrsh && encoderbefore > 1
            encoderjump = encodercurrent;
            break;
        end
    end
end
s1 = data1(encoderjump:encoderjump + encoderXcorrWindow,1);

for i=2:length(test)
    s2 = test{i}(encoderjump: encoderjump + encoderXcorrWindow,1);
    [r,lag] = xcorr(s1,s2);
    [~,I] = max(abs(r));
    encoderjumpavg(i-1) = round(-lag(I));
end

% --- Executes on button press in Align_encoder.
function Align_encoder_Callback(hObject, eventdata, handles)
% hObject    handle to Align_encoder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
accelFile = get(handles.Import_accel,'UserData');
if isempty(accelFile); return; end
AccelFile = accelFile.AccelFile;
AccelPath = accelFile.AccelPath;
AccelFreq = str2double(get(handles.Freq_textbox,'String'));
if isnan(AccelFreq)
    error('Frequency invalid.');
end

cardSN = [1 2 3 4];
chCnt = 4;
encoderVelDis = 3;
[~,fnameraw,fnameext] = fileparts(AccelFile);

test = cell(1,length(cardSN));
try
    for i = 1:length(cardSN)
        test{i} = dlmread([AccelPath, fnameraw(1:end-1),...
            num2str(cardSN(i)), fnameext]);
    end
catch ME
    if i == 1
        try
            test{1} = dlmread([AccelPath,AccelFile]);
            for k = 2:length(cardSN)
                test{k} = test{1};
            end
        catch innerME
            error('Read Accel failed.');
        end
    else
        warning('Not all four cards available.');
        for j = i:length(cardSN)
            test{j} = test{1};
        end
    end
end
set(hObject,'UserData',test);

encoderjumpavg = getCardsTimeDiff(test);

%create time col vector
Acceltime_s = 1/AccelFreq * (1:length(test{1}))';

%create accel data mat
Accelrawdatamat = zeros(length(test{1}),...
    length(cardSN)*(chCnt+1)+encoderVelDis);
Accelrawdatamat(:,1:(chCnt+1)) = [Acceltime_s test{1}];

for i=2:length(cardSN)
    Accelrawdatamat(:,(i-1)*(chCnt+1)+1:i*(chCnt+1)) = ...
        [Acceltime_s-encoderjumpavg(i-1)/AccelFreq test{i}]; % apply time shift
end

%calc encodercounts and velocity

encoderstep = str2double(get(handles.EncoderStep,'String'));
intvl = str2double(get(handles.VelocityWindow,'String'));
if isnan(encoderstep) || isnan(intvl)
    error('Encoder step or velocity window invalid.');
end

[counter,d,v] = getCounterFromEncoder(test{1},encoderstep,AccelFreq,intvl);

Accelrawdatamat(:,end-encoderVelDis+1:end)...
    = [counter d v];
set(handles.Calc_vel,'UserData',Accelrawdatamat);
guidata(hObject,handles);

% align encoder
if isempty(get(handles.Import_c40,'UserData'))
    set(handles.ReAlign,'Enable','off');
    set(handles.Fine_tune_encoder,'Enable','off');
    set(handles.Import_c40,'Enable','off');
    msgbox('Import accel successful. Note there is no c40 file.');
    return;
end
ReAlign_Callback(handles.ReAlign, eventdata, handles);



function [counter,d,v] = getCounterFromEncoder(test,encoderstep,freq,intvl)
lenPerSector = encoderstep; % in meters
f = freq;    % frequency in Hz

count = 0;
encoderShift = 3;
temp = zeros(1,length(test));
shiftedEncoder = test(:,1) - encoderShift;
for j=1:length(test(:,1))-1
    if (shiftedEncoder(j))*(shiftedEncoder(j+1))<0
        count = count + 1;
    end
    temp(j+1) = count;
end

% averaging method of calc velocity
tempv = zeros(1,length(test));
denom = intvl/(f*lenPerSector);
for k=intvl+1:length(temp)-intvl
    tempv(1,k) = -(temp(k-intvl)-temp(k+intvl))/denom;
end

d = temp' * lenPerSector;
counter = temp';
v = tempv';


% --- Executes on slider movement.
function Fine_tune_encoder_Callback(hObject, eventdata, handles)
% hObject    handle to Fine_tune_encoder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Accelrawdatamat = get(handles.Calc_vel,'UserData');
if isempty(Accelrawdatamat); return; end
% offset time vector
offset = get(hObject,'Value');
offset1 = get(hObject,'UserData');  % offset1 is the prev offset
set(hObject,'UserData',offset);
trueoffset = offset-offset1;
timeshift = get(handles.ReAlign,'UserData');
timeshift = timeshift + trueoffset;
set(handles.ReAlign,'UserData',timeshift);

for i=1:4
    Accelrawdatamat(:,5*i-4)=Accelrawdatamat(:,5*i-4)+trueoffset;
end
set(handles.Calc_vel,'UserData',Accelrawdatamat);

% plot new encoder line
lowplot = get(handles.axes2,'UserData');
if isempty(lowplot); return; end
if lowplot ~= 0
    delete(lowplot);
end
hold(handles.axes2,'on');
lowplot = plot(handles.axes2,Accelrawdatamat(:,1),...
    Accelrawdatamat(:,2),'r');
hold(handles.axes2,'off');

set(handles.axes2,'UserData',lowplot);
guidata(hObject,handles);
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function Fine_tune_encoder_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Fine_tune_encoder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on selection change in ydatalist.
function ydatalist_Callback(hObject, eventdata, handles)
% hObject    handle to ydatalist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
str = cellstr(get(hObject,'String'));
val = get(hObject,'Value');
% Set current data to the selected data set.
switch str{val}
    case 'Distance_m'
        ylistsn = 22;
    case 'Velocity'
        ylistsn = 23;
    case 'Normal Stress'
        ylistsn = 24;
    case 'Shear Stress'
        ylistsn = 25;
    case 'Dilation'
        ylistsn = 26;
    case 'Friction'
        ylistsn = 27;
    case 'Temperature'
        ylistsn = 28;
    case 'EncoderCounts'
        ylistsn = 29;
    case 'EncoderVolts'
        ylistsn = 30;
    case 'eddy1_ROGA'
        ylistsn = 31;
    case 'eddy2_ROGA'
        ylistsn = 32;
    case 'eddy3_ROGA'
        ylistsn = 33;
    case 'eddy4_ROGA'
        ylistsn = 34;
    case 'MER'
        ylistsn = 40;
    case 'MotorSpeed'
        ylistsn = 2;
    case 'MotorTorque'
        ylistsn = 14;
    case 'RequestedControlVolts'
        ylistsn = 19;
    case 'RequestedProgramVolts'
        ylistsn = 20;
    case 'Card1 Ch2 XX'
        ylistsn = 103;
    case 'Card1 Ch4 YY'
        ylistsn = 105;
    case 'Card2 Ch2 XX'
        ylistsn = 108;
    case 'Card2 Ch4 YY'
        ylistsn = 110;
    case 'Card3 Ch2 XX'
        ylistsn = 113;
    case 'Card3 Ch4 YY'
        ylistsn = 115;
    case 'Card4 Ch2 XX'
        ylistsn = 118;
    case 'Card4 Ch4 YY'
        ylistsn = 120;
    case 'Card1 Ch3 XY'
        ylistsn = 104;
    case 'Card2 Ch3 XY'
        ylistsn = 109;
    case 'Card3 Ch3 XY'
        ylistsn = 114;
    case 'Card4 Ch3 XY'
        ylistsn = 119;
    case 'Calculated EncoderCounts'
        ylistsn = 121;
    case 'Calculated Distance'
        ylistsn = 122;
    case 'Calculated Velocity'
        ylistsn = 123;
    case 'Encoder_Volts'
        ylistsn = 102;
    case 'None'
        ylistsn = 0;
        
end
set(hObject,'UserData',ylistsn);
% Save the handles structure.
guidata(hObject,handles);
% Hints: contents = cellstr(get(hObject,'String')) returns ydatalist contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ydatalist


% --- Executes during object creation, after setting all properties.
function ydatalist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ydatalist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject, 'String',...
    {'Distance_m';'Velocity';'Normal Stress';'Shear Stress';'Dilation';'Friction';'Temperature';'EncoderCounts';'EncoderVolts';'eddy1_ROGA';'eddy2_ROGA';'eddy3_ROGA';'eddy4_ROGA';'MER';'MotorSpeed';'MotorTorque';'RequestedControlVolts';'RequestedProgramVolts';'Card1 Ch2 XX';'Card1 Ch3 XY';'Card1 Ch4 YY';'Card2 Ch2 XX';'Card2 Ch3 XY';'Card2 Ch4 YY';'Card3 Ch2 XX';'Card3 Ch3 XY';'Card3 Ch4 YY';'Card4 Ch2 XX';'Card4 Ch3 XY';'Card4 Ch4 YY';'Calculated EncoderCounts';'Calculated Distance';'Calculated Velocity';'Encoder_Volts';'None'});
set(hObject, 'Value', 25);
guidata(hObject,handles);


% --- Executes on selection change in plotmethod.
function plotmethod_Callback(hObject, eventdata, handles)
% hObject    handle to plotmethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
str = cellstr(get(hObject,'String'));
val = get(hObject,'Value');
switch str{val}
    case 'linear'
        set(handles.axes1,'XScale','linear');
        set(handles.axes1,'YScale','linear');
    case 'loglog'
        set(handles.axes1,'XScale','log');
        set(handles.axes1,'YScale','log');
    case 'loglinear'
        set(handles.axes1,'XScale','log');
        set(handles.axes1,'YScale','linear');
    case 'linearlog'
        set(handles.axes1,'XScale','linear');
        set(handles.axes1,'YScale','log');
end
guidata(hObject,handles);

% Hints: contents = cellstr(get(hObject,'String')) returns plotmethod contents as cell array
%        contents{get(hObject,'Value')} returns selected item from plotmethod


% --- Executes during object creation, after setting all properties.
function plotmethod_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plotmethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in eventlist.
function eventlist_Callback(hObject, eventdata, handles)
% hObject    handle to eventlist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
T = get(hObject,'UserData');
if isempty(T); 
    set(hObject,'String','Event list');
    set(hObject,'Value',1);
    return;
end
strings = cellstr(num2str(T.SN));
set(hObject,'String',strings);
if get(hObject,'Value') > size(T,1)
    set(hObject,'Value',size(T,1));
end
AccelCh_Callback(handles.AccelCh, eventdata, handles);
% Hints: contents = cellstr(get(hObject,'String')) returns eventlist contents as cell array
%        contents{get(hObject,'Value')} returns selected item from eventlist


% --- Executes during object creation, after setting all properties.
function eventlist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to eventlist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

eventsDB = table();
set(hObject,'UserData',eventsDB);
set(hObject,'String','Event list');
set(hObject,'Value',1);
guidata(hObject,handles);


% --- Executes on button press in add_event.
function add_event_Callback(hObject, eventdata, handles)
% hObject    handle to add_event (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
xl = get(handles.axes1,'XLim');
picktime = clock;
SN = getSN(handles);
Arrival = zeros(3,4);
cellevent = {'SN','pickTime','startEnd','Arrivals';...
    SN,picktime,xl,Arrival};
Tevent = cell2table(cellevent(2:end,:));
Tevent.Properties.VariableNames = cellevent(1,:);
T = get(handles.eventlist,'UserData');
Tnew = [T;Tevent];
set(handles.eventlist,'UserData',Tnew);
set(handles.eventlist,'Value',length(Tnew.SN));
eventlist_Callback(handles.eventlist, eventdata, handles);
guidata(hObject,handles);

function SN = getSN(handles)
T = get(handles.eventlist,'UserData');
if isempty(T)
    SN = 1;
else
    SN = max(T.SN) + 1;
end

function picks = getcursor(handles)
try
    dcm_obj = datacursormode(handles.figure1);
    c_info = getCursorInfo(dcm_obj);
    picks = c_info.Position(1);
catch ME
    msgbox('Please repick.');
    error('Please repick.');
end


% --- Executes on button press in Delete_event.
function Delete_event_Callback(hObject, eventdata, handles)
% hObject    handle to Delete_event (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
T = get(handles.eventlist,'UserData');
if isempty(T); return; end
eventNo = get(handles.eventlist,'Value');
T(eventNo,:) = [];
set(handles.eventlist,'UserData',T);
guidata(hObject,handles);
eventlist_Callback(handles.eventlist, eventdata, handles);

% --- Executes on button press in Spectral.
function Spectral_Callback(hObject, eventdata, handles)
% hObject    handle to Spectral (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
plotSpec(handles,0);

function [f,P1] = plotSpec(handles,withEGF)
T = get(handles.eventlist,'UserData');
if isempty(T); return; end
Accel = get(handles.Calc_vel,'UserData');
if isempty(Accel); return; end
eventNo = get(handles.eventlist,'Value');
i = get(handles.AccelCh,'Value');
card = floor((i-1)/3)+1;
ch = mod(i-1,3)+1;
Fs = str2double(get(handles.Freq_textbox,'String'));
if isnan(Fs); return; end
SNs = [
     3     8    13    18
     4     9    14    19
     5    10    15    20
     ];
ind1 = find(Accel(:,SNs(1,card) - 2)>=T.startEnd(eventNo,1),1);
ind2 = find(Accel(:,SNs(1,card) - 2)>=T.startEnd(eventNo,2),1);
% t = Accel(ind1:ind2,SNs(1,card) - 2);
sig = Accel(ind1:ind2,SNs(ch,card));
if withEGF
    EGFSN = get(handles.Assign_EGF,'UserData');
    if isempty(EGFSN); return; end
    EGFNo = find(T.SN == EGFSN,1);
    if isempty(EGFNo)
        error('EGF does not exist.');
    end
    ind1 = find(Accel(:,SNs(1,card) - 2)>=T.startEnd(EGFNo,1),1);
    ind2 = find(Accel(:,SNs(1,card) - 2)>=T.startEnd(EGFNo,2),1);
    EGFSig = Accel(ind1:ind2,SNs(ch,card));
    N = max(length(sig),length(EGFSig));
    [f,P1Sig] = calcFFT(sig,Fs,N);
    [~,P1EGFSig] = calcFFT(EGFSig,Fs,N);
    loglog(handles.axes2,f,P1Sig./P1EGFSig);
else
    N = length(sig);
    [f,P1] = calcFFT(sig,Fs,N);
    loglog(handles.axes2,f,P1);
end


function [f,P1] = calcFFT(sig,Fs,N)
x = sig-mean(sig);
y = fft(x,N);
P2 = abs(y/N);
P1 = P2(1:floor(N/2)+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:floor(N/2))/N;

% --- Executes on button press in Assign_EGF.
function Assign_EGF_Callback(hObject, eventdata, handles)
% hObject    handle to Assign_EGF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~iscell(get(handles.eventlist,'String'))
    return;
end
EGFNo = get(handles.eventlist,'Value');
strings = get(handles.eventlist,'String');
SN = strings{EGFNo};
SNNo = str2double(SN);
set(hObject,'UserData',SNNo);
set(hObject,'String',['EGF: ',SN]);
guidata(hObject,handles);

% --- Executes on button press in Spec_ratio.
function Spec_ratio_Callback(hObject, eventdata, handles)
% hObject    handle to Spec_ratio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
plotSpec(handles,1);

% --- Executes on button press in Debug.
function Debug_Callback(hObject, eventdata, handles)
% hObject    handle to Debug (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
keyboard;


% --- Executes on button press in Zoom.
function Zoom_Callback(hObject, eventdata, handles)
% hObject    handle to Zoom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
x1=1e-6;
set(hObject,'UserData',[]);
while (x1 > 0)
    [x,~,ibut] = ginput (2);
    
    if (ibut == 1)
        xlim([min(x) max(x)]);
    end
    
    if (ibut == 2)
        xlim([mean(x)-2*abs(x(2)-x(1)) mean(x)+10*abs(x(2)-x(1))]);
    end
    
    if (ibut(1) == 1 && ibut(2) == 3)
        xlim('auto');
        ylim('auto');
    end
    
    if (ibut == 3)
        set(hObject,'UserData',x);
        set(hObject,'String','Picked!');
    end
    
    if (ibut(1) == 3 && ibut(2) == 1)
        x1=0;
    end
end
set(hObject,'String','Zoom');


% --- Executes on button press in Full_res.
function Full_res_Callback(hObject, eventdata, handles)
% hObject    handle to Full_res (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Plot_Callback(handles.Plot, eventdata, handles);
% Hint: get(hObject,'Value') returns toggle state of Full_res


% --- Executes on button press in save.
function save_Callback(hObject, eventdata, handles)
% hObject    handle to save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Accel = get(handles.Calc_vel,'UserData');
ROGA = get(handles.Plot,'UserData');
eventsDB = get(handles.eventlist,'UserData');
test = get(handles.Align_encoder,'UserData');
AccelFile = get(handles.Import_accel,'UserData');
ROGAFile = get(handles.Import_c40,'UserData');
try
    [~,dryname,~] = fileparts(AccelFile.AccelFile);
    save([AccelFile.AccelPath, 'DB_', dryname(1:end-2), '.mat'],...
        'Accel','ROGA','eventsDB','test','AccelFile','ROGAFile');
    msgbox('Save successful!');
catch ME
    msgbox('Save failed.');
end



% --- Executes on button press in Load.
function Load_Callback(hObject, eventdata, handles)
% hObject    handle to Load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[file,path] = uigetfile(...
    'C:\Users\User\Documents\Test files results on fast cards\*.mat');
if file == 0; return; end
clearUserData(handles);
toInput = {'Accel','ROGA','eventsDB','test','AccelFile','ROGAFile'};
s = load([path,file],toInput{:});
try
    set(handles.Calc_vel,'UserData',s.Accel);
    set(handles.eventlist,'UserData',s.eventsDB);
    set(handles.Plot,'UserData',s.ROGA);
    set(handles.Align_encoder,'UserData',s.test,'Enable','off');
    AccelPathFile = s.AccelFile;
    set(handles.Import_accel,'UserData',AccelPathFile);
    set(handles.Import_accel,'TooltipString',[AccelPathFile.AccelPath,AccelPathFile.AccelFile]);
    set(handles.Import_c40,'UserData',s.ROGAFile,'Enable','off');
    set(handles.Fine_tune_encoder,'Enable','off');
    set(handles.ReAlign,'Enable','off');
    guidata(hObject,handles);
    eventlist_Callback(handles.eventlist, eventdata, handles);
    msgbox('Load successful!');
catch ME
    msgbox('Load failed.');
end


% --- Executes on selection change in ydatalist2.
function ydatalist2_Callback(hObject, eventdata, handles)
% hObject    handle to ydatalist2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
str = cellstr(get(hObject,'String'));
val = get(hObject,'Value');
% Set current data to the selected data set.
switch str{val}
    case 'Distance_m'
        ylistsn = 22;
    case 'Velocity'
        ylistsn = 23;
    case 'Normal Stress'
        ylistsn = 24;
    case 'Shear Stress'
        ylistsn = 25;
    case 'Dilation'
        ylistsn = 26;
    case 'Friction'
        ylistsn = 27;
    case 'Temperature'
        ylistsn = 28;
    case 'EncoderCounts'
        ylistsn = 29;
    case 'EncoderVolts'
        ylistsn = 30;
    case 'eddy1_ROGA'
        ylistsn = 31;
    case 'eddy2_ROGA'
        ylistsn = 32;
    case 'eddy3_ROGA'
        ylistsn = 33;
    case 'eddy4_ROGA'
        ylistsn = 34;
    case 'MER'
        ylistsn = 40;
    case 'MotorSpeed'
        ylistsn = 2;
    case 'MotorTorque'
        ylistsn = 14;
    case 'RequestedControlVolts'
        ylistsn = 19;
    case 'RequestedProgramVolts'
        ylistsn = 20;
    case 'Card1 Ch2 XX'
        ylistsn = 103;
    case 'Card1 Ch4 YY'
        ylistsn = 105;
    case 'Card2 Ch2 XX'
        ylistsn = 108;
    case 'Card2 Ch4 YY'
        ylistsn = 110;
    case 'Card3 Ch2 XX'
        ylistsn = 113;
    case 'Card3 Ch4 YY'
        ylistsn = 115;
    case 'Card4 Ch2 XX'
        ylistsn = 118;
    case 'Card4 Ch4 YY'
        ylistsn = 120;
    case 'Card1 Ch3 XY'
        ylistsn = 104;
    case 'Card2 Ch3 XY'
        ylistsn = 109;
    case 'Card3 Ch3 XY'
        ylistsn = 114;
    case 'Card4 Ch3 XY'
        ylistsn = 119;
    case 'Calculated EncoderCounts'
        ylistsn = 121;
    case 'Calculated Distance'
        ylistsn = 122;
    case 'Calculated Velocity'
        ylistsn = 123;
    case 'Encoder_Volts'
        ylistsn = 102;
    case 'None'
        ylistsn = 0;
        
end
set(hObject,'UserData',ylistsn);
% Save the handles structure.
guidata(hObject,handles);
% Hints: contents = cellstr(get(hObject,'String')) returns ydatalist2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ydatalist2


% --- Executes during object creation, after setting all properties.
function ydatalist2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ydatalist2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject, 'String',...
    {'Distance_m';'Velocity';'Normal Stress';'Shear Stress';'Dilation';'Friction';'Temperature';'EncoderCounts';'EncoderVolts';'eddy1_ROGA';'eddy2_ROGA';'eddy3_ROGA';'eddy4_ROGA';'MER';'MotorSpeed';'MotorTorque';'RequestedControlVolts';'RequestedProgramVolts';'Card1 Ch2 XX';'Card1 Ch3 XY';'Card1 Ch4 YY';'Card2 Ch2 XX';'Card2 Ch3 XY';'Card2 Ch4 YY';'Card3 Ch2 XX';'Card3 Ch3 XY';'Card3 Ch4 YY';'Card4 Ch2 XX';'Card4 Ch3 XY';'Card4 Ch4 YY';'Calculated EncoderCounts';'Calculated Distance';'Calculated Velocity';'Encoder_Volts';'None'});
set(hObject, 'Value', 6);
guidata(hObject,handles);


% --- Executes on selection change in ydatalist3.
function ydatalist3_Callback(hObject, eventdata, handles)
% hObject    handle to ydatalist3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
str = cellstr(get(hObject,'String'));
val = get(hObject,'Value');
% Set current data to the selected data set.
switch str{val}
    case 'Distance_m'
        ylistsn = 22;
    case 'Velocity'
        ylistsn = 23;
    case 'Normal Stress'
        ylistsn = 24;
    case 'Shear Stress'
        ylistsn = 25;
    case 'Dilation'
        ylistsn = 26;
    case 'Friction'
        ylistsn = 27;
    case 'Temperature'
        ylistsn = 28;
    case 'EncoderCounts'
        ylistsn = 29;
    case 'EncoderVolts'
        ylistsn = 30;
    case 'eddy1_ROGA'
        ylistsn = 31;
    case 'eddy2_ROGA'
        ylistsn = 32;
    case 'eddy3_ROGA'
        ylistsn = 33;
    case 'eddy4_ROGA'
        ylistsn = 34;
    case 'MER'
        ylistsn = 40;
    case 'MotorSpeed'
        ylistsn = 2;
    case 'MotorTorque'
        ylistsn = 14;
    case 'RequestedControlVolts'
        ylistsn = 19;
    case 'RequestedProgramVolts'
        ylistsn = 20;
    case 'Card1 Ch2 XX'
        ylistsn = 103;
    case 'Card1 Ch4 YY'
        ylistsn = 105;
    case 'Card2 Ch2 XX'
        ylistsn = 108;
    case 'Card2 Ch4 YY'
        ylistsn = 110;
    case 'Card3 Ch2 XX'
        ylistsn = 113;
    case 'Card3 Ch4 YY'
        ylistsn = 115;
    case 'Card4 Ch2 XX'
        ylistsn = 118;
    case 'Card4 Ch4 YY'
        ylistsn = 120;
    case 'Card1 Ch3 XY'
        ylistsn = 104;
    case 'Card2 Ch3 XY'
        ylistsn = 109;
    case 'Card3 Ch3 XY'
        ylistsn = 114;
    case 'Card4 Ch3 XY'
        ylistsn = 119;
    case 'Calculated EncoderCounts'
        ylistsn = 121;
    case 'Calculated Distance'
        ylistsn = 122;
    case 'Calculated Velocity'
        ylistsn = 123;
    case 'Encoder_Volts'
        ylistsn = 102;
    case 'None'
        ylistsn = 0;
        
end
set(hObject,'UserData',ylistsn);
% Save the handles structure.
guidata(hObject,handles);
% Hints: contents = cellstr(get(hObject,'String')) returns ydatalist3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ydatalist3


% --- Executes during object creation, after setting all properties.
function ydatalist3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ydatalist3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject, 'String',...
    {'Distance_m';'Velocity';'Normal Stress';'Shear Stress';'Dilation';'Friction';'Temperature';'EncoderCounts';'EncoderVolts';'eddy1_ROGA';'eddy2_ROGA';'eddy3_ROGA';'eddy4_ROGA';'MER';'MotorSpeed';'MotorTorque';'RequestedControlVolts';'RequestedProgramVolts';'Card1 Ch2 XX';'Card1 Ch3 XY';'Card1 Ch4 YY';'Card2 Ch2 XX';'Card2 Ch3 XY';'Card2 Ch4 YY';'Card3 Ch2 XX';'Card3 Ch3 XY';'Card3 Ch4 YY';'Card4 Ch2 XX';'Card4 Ch3 XY';'Card4 Ch4 YY';'Calculated EncoderCounts';'Calculated Distance';'Calculated Velocity';'Encoder_Volts';'None'});
set(hObject, 'Value', 26);
guidata(hObject,handles);

% --- Executes on selection change in ydatalist4.
function ydatalist4_Callback(hObject, eventdata, handles)
% hObject    handle to ydatalist4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
str = cellstr(get(hObject,'String'));
val = get(hObject,'Value');
% Set current data to the selected data set.
switch str{val}
    case 'Distance_m'
        ylistsn = 22;
    case 'Velocity'
        ylistsn = 23;
    case 'Normal Stress'
        ylistsn = 24;
    case 'Shear Stress'
        ylistsn = 25;
    case 'Dilation'
        ylistsn = 26;
    case 'Friction'
        ylistsn = 27;
    case 'Temperature'
        ylistsn = 28;
    case 'EncoderCounts'
        ylistsn = 29;
    case 'EncoderVolts'
        ylistsn = 30;
    case 'eddy1_ROGA'
        ylistsn = 31;
    case 'eddy2_ROGA'
        ylistsn = 32;
    case 'eddy3_ROGA'
        ylistsn = 33;
    case 'eddy4_ROGA'
        ylistsn = 34;
    case 'MER'
        ylistsn = 40;
    case 'MotorSpeed'
        ylistsn = 2;
    case 'MotorTorque'
        ylistsn = 14;
    case 'RequestedControlVolts'
        ylistsn = 19;
    case 'RequestedProgramVolts'
        ylistsn = 20;
    case 'Card1 Ch2 XX'
        ylistsn = 103;
    case 'Card1 Ch4 YY'
        ylistsn = 105;
    case 'Card2 Ch2 XX'
        ylistsn = 108;
    case 'Card2 Ch4 YY'
        ylistsn = 110;
    case 'Card3 Ch2 XX'
        ylistsn = 113;
    case 'Card3 Ch4 YY'
        ylistsn = 115;
    case 'Card4 Ch2 XX'
        ylistsn = 118;
    case 'Card4 Ch4 YY'
        ylistsn = 120;
    case 'Card1 Ch3 XY'
        ylistsn = 104;
    case 'Card2 Ch3 XY'
        ylistsn = 109;
    case 'Card3 Ch3 XY'
        ylistsn = 114;
    case 'Card4 Ch3 XY'
        ylistsn = 119;
    case 'Calculated EncoderCounts'
        ylistsn = 121;
    case 'Calculated Distance'
        ylistsn = 122;
    case 'Calculated Velocity'
        ylistsn = 123;
    case 'Encoder_Volts'
        ylistsn = 102;
    case 'None'
        ylistsn = 0;

    end
set(hObject,'UserData',ylistsn);
% Save the handles structure.
guidata(hObject,handles);
% Hints: contents = cellstr(get(hObject,'String')) returns ydatalist4 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ydatalist4


% --- Executes during object creation, after setting all properties.
function ydatalist4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ydatalist4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject, 'String',...
    {'Distance_m';'Velocity';'Normal Stress';'Shear Stress';'Dilation';'Friction';'Temperature';'EncoderCounts';'EncoderVolts';'eddy1_ROGA';'eddy2_ROGA';'eddy3_ROGA';'eddy4_ROGA';'MER';'MotorSpeed';'MotorTorque';'RequestedControlVolts';'RequestedProgramVolts';'Card1 Ch2 XX';'Card1 Ch3 XY';'Card1 Ch4 YY';'Card2 Ch2 XX';'Card2 Ch3 XY';'Card2 Ch4 YY';'Card3 Ch2 XX';'Card3 Ch3 XY';'Card3 Ch4 YY';'Card4 Ch2 XX';'Card4 Ch3 XY';'Card4 Ch4 YY';'Calculated EncoderCounts';'Calculated Distance';'Calculated Velocity';'Encoder_Volts';'None'});
set(hObject, 'Value', 27);
guidata(hObject,handles);



function EncoderStep_Callback(hObject, eventdata, handles)
% hObject    handle to EncoderStep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EncoderStep as text
%        str2double(get(hObject,'String')) returns contents of EncoderStep as a double


% --- Executes during object creation, after setting all properties.
function EncoderStep_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EncoderStep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function VelocityWindow_Callback(hObject, eventdata, handles)
% hObject    handle to VelocityWindow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of VelocityWindow as text
%        str2double(get(hObject,'String')) returns contents of VelocityWindow as a double


% --- Executes during object creation, after setting all properties.
function VelocityWindow_CreateFcn(hObject, eventdata, handles)
% hObject    handle to VelocityWindow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Calc_vel.
function Calc_vel_Callback(hObject, eventdata, handles)
% hObject    handle to Calc_vel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
intvl = str2double(get(handles.VelocityWindow,'String'));
if isnan(intvl); error('Invalid Vel Window.'); end
Accelrawdatamat = get(hObject,'UserData');
if isempty(Accelrawdatamat); return; end
encoderstep = str2double(get(handles.EncoderStep,'String'));
if isnan(encoderstep); error('Invalid encoder step.'); end
test = get(handles.Align_encoder,'UserData');
if isempty(test); warning('no source file.'); return; end
freq = str2double(get(handles.Freq_textbox,'String'));
if isnan(freq); error('Invalid Freq.'); end

encoderVelDis = 3;
[counter,d,v] = getCounterFromEncoder(test{1},encoderstep,freq,intvl);
Accelrawdatamat(:,end-encoderVelDis+1:end)...
    = [counter d v];
set(handles.Calc_vel,'UserData',Accelrawdatamat);

guidata(hObject,handles);
Plot_Callback(handles.Plot, eventdata, handles);


% --- Executes on button press in Scaled.
function Scaled_Callback(hObject, eventdata, handles)
% hObject    handle to Scaled (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Plot_Callback(handles.Plot, eventdata, handles);
% Hint: get(hObject,'Value') returns toggle state of Scaled


% --- Executes on button press in Plot_NSE.
function Plot_NSE_Callback(hObject, eventdata, handles)
% hObject    handle to Plot_NSE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Plot_Callback(handles.Plot, eventdata, handles);
T = get(handles.eventlist,'UserData');
if isempty(T); return; end

hold(handles.axes1,'on');
yl = get(handles.axes1,'YLim');
SN = T.SN;
start = T.startEnd(:,1);
for i = 1:length(SN)
    plot(handles.axes1,ones(1,2)*start(i),yl,...
        'Color','black');
    text(start(i),getPlotHeight(i,yl),num2str(SN(i)),...
        'FontSize',20);
end
hold(handles.axes1,'off');

function h = getPlotHeight(i,yl)
height = (mod(i,20)+1)/20;
h = (1-height)*yl(1) + height*yl(2);


% --- Executes on button press in ReAlign.
function ReAlign_Callback(hObject, eventdata, handles)
% hObject    handle to ReAlign (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ROGAc40mat = get(handles.Plot,'UserData');
Accelrawdatamat = get(handles.Calc_vel,'UserData');
if isempty(ROGAc40mat) || isempty(Accelrawdatamat)
    return;
end

%plot encoder on both plot and align
cla(handles.axes1);
hold(handles.axes1,'on');
plot(handles.axes1,Accelrawdatamat(:,1),Accelrawdatamat(:,2)); % encoder on Accel
plot(handles.axes1,Accelrawdatamat(:,1),Accelrawdatamat(:,end-1)); % encoder d
plot(handles.axes1,Accelrawdatamat(:,1),Accelrawdatamat(:,end)); % encoder v
hold(handles.axes1,'off');

cla(handles.axes2);
hold(handles.axes2,'on');
plot(handles.axes2,ROGAc40mat(:,21),ROGAc40mat(:,13)); % encoder on ROGA
plot(handles.axes2,ROGAc40mat(:,21),ROGAc40mat(:,22)); % d on ROGA
plot(handles.axes2,ROGAc40mat(:,21),ROGAc40mat(:,23)); % velocity on ROGA
hold(handles.axes2,'off');

% mouse click to select starting point on high and low freq plot to
% align encoder

Zoom_Callback(handles.Zoom, eventdata, handles);
x = get(handles.Zoom,'UserData');
if isempty(x); error('No picking performed!'); end
AccelFreq = str2double(get(handles.Freq_textbox,'String'));
if isnan(AccelFreq); error('Invalid freq.'); end
ROGAFreq = get(handles.Freq_textbox,'UserData');

encoderlow = ROGAc40mat(:,13);
event = Accelrawdatamat(:,2);

lowmark = x(1); % first right click on low freq
highmark = x(2);    % second right click on high freq

jlow=0;
kk=round(lowmark*ROGAFreq);
while jlow==0
    if(abs(round(encoderlow(kk))-round(encoderlow(kk+1)))<2)
        kk=kk+1;
    else jlow=kk-1;
    end
    if kk==length(encoderlow);
        break;
    end
end

kk=round((highmark-Accelrawdatamat(1,1))*AccelFreq);
jhigh=0;
while jhigh==0
    if(abs(round(event(kk,1))-round(event(kk+1,1)))<2)
        kk=kk+1;
    else jhigh=kk-1;
    end
end

timeshift = (jlow-1)/ROGAFreq + ROGAc40mat(1,21)...
    -(jhigh-1)/AccelFreq-Accelrawdatamat(1,1);
set(hObject,'UserData',timeshift);
for i = 1:4
    Accelrawdatamat(:,5*i-4)=Accelrawdatamat(:,5*i-4)+timeshift;
end
set(handles.Calc_vel,'UserData',Accelrawdatamat);
hold(handles.axes2,'on');
lowplot = plot(handles.axes2,Accelrawdatamat(:,1),...
    event(:,1),'r');
set(handles.axes2,'UserData',lowplot);
hold(handles.axes2,'off');

guidata(hObject,handles);


% --- Executes when figure1 is resized.
function figure1_SizeChangedFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function Full_res_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Full_res (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
set(hObject,'Value',get(hObject,'Min'));


% --- Executes during object creation, after setting all properties.
function add_event_CreateFcn(hObject, eventdata, handles)
% hObject    handle to add_event (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on selection change in AccelCh.
function AccelCh_Callback(hObject, eventdata, handles)
% hObject    handle to AccelCh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
eventNo = get(handles.eventlist,'Value');
i = get(hObject,'Value');
card = floor((i-1)/3)+1;
ch = mod(i-1,3)+1;
T = get(handles.eventlist,'UserData');
if isempty(T); return; end
lastItem = 35;
set(handles.ydatalist2,'Value',lastItem);
set(handles.ydatalist3,'Value',lastItem);
set(handles.ydatalist4,'Value',lastItem);
% set(handles.Full_res,'Value',0);
% set(handles.Scaled,'Value',1);
SNs = [19 21 23 25;
       27 28 29 30;
       20 22 24 26];
set(handles.ydatalist,'Value',SNs(ch,card));
Plot_Callback(handles.Plot, eventdata, handles);
yl = get(handles.axes1,'YLim');
% xl = get(handles.axes1,'XLim');
set(handles.axes1,'XLim',T.startEnd(eventNo,:));
zoom reset;

% update unpicked channels
unpicked = find(0==cell2mat(T.Arrivals(eventNo)));
toPrint = cell(1,length(unpicked)+1);
toPrint{1,1} = 'Unpicked Channels:';
for i = 1:length(unpicked)
    ich = mod(unpicked(i)-1,3)+1;
    icard = floor((unpicked(i)-1)/3)+1;
    toPrint{1,i+1} = ['Accel ',num2str(icard),...
        ' Ch ', num2str(ich)];
end
set(handles.Unpicked,'String',toPrint);

% plot a line to show the pick
pick = T.Arrivals{eventNo}(ch,card);
if pick ~= 0
    hold(handles.axes1,'on');
    plot(handles.axes1,ones(1,2)*pick,yl,'r');
    % plot(handles.axes1,xl,ones(1,2)*pick,'r');
    hold(handles.axes1,'off');
end

guidata(hObject,handles);
Spectral_Callback(handles.Spectral, eventdata, handles);

% Hints: contents = cellstr(get(hObject,'String')) returns AccelCh contents as cell array
%        contents{get(hObject,'Value')} returns selected item from AccelCh


% --- Executes during object creation, after setting all properties.
function AccelCh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to AccelCh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
liststrings = cell(1,3*4);
for i = 1:3*4
    liststrings{i} = ['Accel ',num2str(floor((i-1)/3)+1),' Ch ',...
        num2str(mod(i-1,3)+1)];
end
set(hObject,'String',liststrings);


% --- Executes on button press in RegPick.
function RegPick_Callback(hObject, eventdata, handles)
% hObject    handle to RegPick (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
T = get(handles.eventlist,'UserData');
if isempty(T); return; end
pick = getcursor(handles);
i = get(handles.AccelCh,'Value');
card = floor((i-1)/3)+1;
ch = mod(i-1,3)+1;
eventNo = get(handles.eventlist,'Value');
T.Arrivals{eventNo}(ch,card) = pick;
set(handles.eventlist,'UserData',T);
AccelCh_Callback(handles.AccelCh, eventdata, handles);
guidata(hObject,handles);


% --- Executes on button press in savetemp.
function savetemp_Callback(hObject, eventdata, handles)
% hObject    handle to savetemp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Accel = get(handles.Calc_vel,'UserData');
ROGA = get(handles.Plot,'UserData');
eventsDB = get(handles.eventlist,'UserData');
if isempty(Accel); return; end

assignin('base','Accel',Accel);
assignin('base','ROGA',ROGA);
assignin('base','eventsDB',eventsDB);


% --- Executes on button press in loadtemp.
function loadtemp_Callback(hObject, eventdata, handles)
% hObject    handle to loadtemp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
    set(handles.Calc_vel,'UserData',evalin('base','Accel'));
    set(handles.Plot,'UserData',evalin('base','ROGA'));
    set(handles.eventlist,'UserData',evalin('base','eventsDB'));
catch ME
    error('load failed.');
end


% --- Executes on button press in chop.
function chop_Callback(hObject, eventdata, handles)
% hObject    handle to chop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Plot_Callback(handles.Plot, eventdata, handles);
% Hint: get(hObject,'Value') returns toggle state of chop


% --- Executes during object creation, after setting all properties.
function Plot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in autoscale.
function autoscale_Callback(hObject, eventdata, handles)
% hObject    handle to autoscale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.axes1,'XLimMode','auto','YLimMode','auto');
zoom reset;


% --- Executes on button press in zoomreset.
function zoomreset_Callback(hObject, eventdata, handles)
% hObject    handle to zoomreset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
zoom reset;

function show = getShowCh(handles)
show = zeros(1,4);
show(1) = get(handles.showCh1,'Value');
show(2) = get(handles.showCh2,'Value');
show(3) = get(handles.showCh3,'Value');
show(4) = get(handles.showCh4,'Value');

% --- Executes on button press in showCh1.
function showCh1_Callback(hObject, eventdata, handles)
% hObject    handle to showCh1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of showCh1


% --- Executes on button press in showCh2.
function showCh2_Callback(hObject, eventdata, handles)
% hObject    handle to showCh2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of showCh2


% --- Executes on button press in showCh3.
function showCh3_Callback(hObject, eventdata, handles)
% hObject    handle to showCh3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of showCh3


% --- Executes on button press in showCh4.
function showCh4_Callback(hObject, eventdata, handles)
% hObject    handle to showCh4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of showCh4


% --- Executes on button press in Smooth.
function Smooth_Callback(hObject, eventdata, handles)
% hObject    handle to Smooth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Plot_Callback(handles.Plot,eventdata,handles);
% Hint: get(hObject,'Value') returns toggle state of Smooth


function offsetAndSlope_Callback(hObject, eventdata, handles)
% hObject    handle to offsetAndSlope (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of offsetAndSlope as text
%        str2double(get(hObject,'String')) returns contents of offsetAndSlope as a double


% --- Executes during object creation, after setting all properties.
function offsetAndSlope_CreateFcn(hObject, eventdata, handles)
% hObject    handle to offsetAndSlope (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
