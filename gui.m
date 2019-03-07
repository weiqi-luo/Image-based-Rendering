function varargout = gui(varargin)
% Last Modified by GUIDE v2.5 11-Sep-2018 17:39:53

%INTRODUCTION: In this gui you can load left image and right image separately 
%and choose p value. You also can choose generate disparity map or just load a
%generated disparity map to save time. Finally you can see the virtual image.
%Clear putton is to reset all the parameter.

%Viel Spa?!

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui_OpeningFcn, ...
                   'gui_OutputFcn',  @gui_OutputFcn, ...
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


% --- Executes just before gui is made visible.
function gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.

% Choose default command line output for gui
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);
axes(handles.axes2);
imshow(imread('img/tum.jpg'))

% --- Outputs from this function are returned to the command line.
function varargout = gui_OutputFcn(hObject, eventdata, handles) 
% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1
str = get(hObject, 'String');
val = get(hObject,'Value');
switch str{val}
    case 'L1'
        I1 = imread('img/L1.jpg');
    case 'L2'
        I1 = imread('img/L2.jpg');
end
handles.popupmenu1 = I1;
assignin('base', 'I1', I1);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2
str = get(hObject, 'String');
val = get(hObject,'Value');
handles.im = val;
switch str{val}
    case 'R1'
        I2 = imread('img/R1.jpg');
    case 'R2'
        I2 = imread('img/R2.jpg');
end
handles.popupmenu2 = I2;
assignin('base', 'I2', I2);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



%Get p value
function edit1_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double
p = str2double(get(hObject,'String'));
if isnan(p)
    p = 0;
    set(hObject,'String','');
    errordlg('Input must be a number', 'Error');
elseif (p>1 || p<0)
    p = 0;
    set(hObject,'String','');
    errordlg('Input must between 0 and 1', 'Error');
end
assignin('base', 'p', p);
handles.edit1 = p;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in radiobutton1. Generate disparity map
function radiobutton1_Callback(hObject, eventdata, handles)
% Hint: get(hObject,'Value') returns toggle state of radiobutton1
disp1 = get(hObject,'Value');
handles.radiobutton1 = disp1;
guidata(hObject, handles);


% --- Executes on button press in radiobutton2.Load disparity map
function radiobutton2_Callback(hObject, eventdata, handles)
% Hint: get(hObject,'Value') returns toggle state of radiobutton2
disp2 = get(hObject,'Value');
handles.radiobutton2 = disp2;
guidata(hObject, handles);




function edit5_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double
ratio = str2double(get(hObject,'String'));
if isnan(ratio)
    ratio = 0;
    set(hObject,'String','');
    errordlg('Input must be a number', 'Error');
elseif (ratio>1 || ratio<0)
    ratio = 0;
    set(hObject,'String','');
    errordlg('Input must between 0 and 1', 'Error');
end
handles.edit5 = ratio;
guidata(hObject, handles);




% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton5. Show virtual image
function pushbutton5_Callback(hObject, eventdata, handles)
p = handles.edit1;
disp2 = handles.radiobutton2;
if disp2 == 1
    load_dis = true;
    ratio = 0.5;
else
    load_dis = false;
    ratio = handles.edit5;
end
I1 = handles.popupmenu1;
I2 = handles.popupmenu2;
if  handles.im == 2
    dispaiy_range = [-500,620];
    choose_img = true;
    Np = 700;
else
    dispaiy_range = [-426,450];
    choose_img = false;
    Np = 1400 ;
end
tic
output_img = free_viewpoint(I1, I2, 'choose_img', choose_img,'load_disparityMap',load_dis, ...
    'p', p, 'down_ratio',ratio ,'disparity_range', dispaiy_range,'Np',Np);
elapsed_time = toc;
time = sprintf('Running time :  %f ',elapsed_time)
set(handles.text3, 'String',time);
guidata(hObject, handles);


% --- Executes on button press in pushbutton6. Clear 
function pushbutton6_Callback(hObject, eventdata, handles)
close(gcbf) 
gui
