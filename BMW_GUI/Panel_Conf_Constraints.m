%% Panel_Conf_Constraints (Configuration: Constraints)
%
% -----------------------
% Programmed by Ma, Tianye
% Under the instruction of Dr. Ku, Yixuan
% Memory, Attention & Cognition (MAC) Lab,
% 11/6/2019
%
% Bug reports or any other feedbacks please contact M.T. (mack_ma2018@outlook.com)
% BMW toolbox: https://github.com/Mack-Ma/Bayesian_Modeling_of_Working_Memory
%

function varargout = Panel_Conf_Constraints(varargin)

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Panel_Conf_Constraints_OpeningFcn, ...
                   'gui_OutputFcn',  @Panel_Conf_Constraints_OutputFcn, ...
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

% --- Executes just before Panel_Conf_Constraints is made visible.
function Panel_Conf_Constraints_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

WindowSize=get(gcf,'Position');
WindowSize(1:2)=[39.8 20];
set(gcf,'Position',WindowSize)

% --- Outputs from this function are returned to the command line.
function varargout = Panel_Conf_Constraints_OutputFcn(hObject, eventdata, handles) 

varargout{1} = handles.output;

% --- Executes on button press in Done_Conf_Constraints.
function Done_Conf_Constraints_Callback(hObject, eventdata, handles)
global Constraints
Constraints.start=eval(['[',get(handles.Blank_Start,'String'),']']);
Constraints.ub=eval(['[',get(handles.Blank_UB,'String'),']']);
Constraints.lb=eval(['[',get(handles.Blank_LB,'String'),']']);
close(gcf);

function Blank_Start_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function Blank_Start_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Blank_UB_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function Blank_UB_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Blank_LB_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function Blank_LB_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
