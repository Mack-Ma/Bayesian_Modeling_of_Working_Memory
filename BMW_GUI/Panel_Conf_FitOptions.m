%% Panel_Conf_FitOptions (Configuration: Fit Options)
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

function varargout = Panel_Conf_FitOptions(varargin)

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Panel_Conf_FitOptions_OpeningFcn, ...
                   'gui_OutputFcn',  @Panel_Conf_FitOptions_OutputFcn, ...
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

% --- Executes just before Panel_Conf_FitOptions is made visible.
function Panel_Conf_FitOptions_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
global AllAlgorithm AllDisplay
% Update handles structure
guidata(hObject, handles);

WindowSize=get(gcf,'Position');
WindowSize(1:2)=[39.8 20];
set(gcf,'Position',WindowSize)

AllAlgorithm={'fmincon: Inferior Point','fmincon: Active Set','fmincon: SQP','fminsearch','CM-AES','BADS','MADS','bayesopt','Genetic Algorithm','Simulated Annealing'};
AllDisplay={'off','iter','notify','final'};
set(handles.List_Algorithm,'String',AllAlgorithm)
set(handles.List_Display,'String',AllDisplay)

% --- Outputs from this function are returned to the command line.
function varargout = Panel_Conf_FitOptions_OutputFcn(hObject, eventdata, handles) 

varargout{1} = handles.output;

% --- Executes on button press in Done_Conf_FitOptions.
function Done_Conf_FitOptions_Callback(hObject, eventdata, handles)
global FitOptions AllAlgorithm AllDisplay
FitOptions.Algorithm=AllAlgorithm{get(handles.List_Algorithm,'Value')};
FitOptions.MaxIter=str2double(get(handles.Blank_MaxIter,'String'));
FitOptions.Display=AllDisplay{get(handles.List_Display,'Value')};
close(gcf);

function List_Algorithm_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function List_Algorithm_CreateFcn(hObject, eventdata, handles)

function Blank_MaxIter_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function Blank_MaxIter_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function List_Display_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function List_Display_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
