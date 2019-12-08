%% Panel_MC (Model Specification)
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

function varargout = Panel_MS(varargin)

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Panel_MS_OpeningFcn, ...
                   'gui_OutputFcn',  @Panel_MS_OutputFcn, ...
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


% --- Executes just before Panel_MS is made visible.
function Panel_MS_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
guidata(hObject, handles);

WindowSize=get(gcf,'Position');
WindowSize(1:2)=[150 25];
set(gcf,'Position',WindowSize)

global MOI
MOI=[];
AllCategory={'','Continuous Recall Models','Change Detection Models','Custom Models'};
set(handles.List_Model,'string',AllCategory); 

% --- Outputs from this function are returned to the command line.
function varargout = Panel_MS_OutputFcn(hObject, eventdata, handles) 

varargout{1} = handles.output;

% --- Executes on selection change in List_Model.
function List_Model_Callback(hObject, eventdata, ~)

if hObject.Value==2
    Panel_MS_CR
elseif hObject.Value==3
    Panel_MS_CD
elseif hObject.Value==4
    Panel_MS_Custom
end

% --- Executes during object creation, after setting all properties.
function List_Model_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in Done_MS.
function Done_MS_Callback(hObject, eventdata, handles)
global MOI
if isempty(MOI)
    error('Model of interest should be specified...')
end
MS_BMW.Model=MOI;
MS_BMW.OutputDir=get(handles.Blank_Output_MS,'String');
MS_BMW.Data=get(handles.Data,'String');
ModelDefinition_BMW(MS_BMW);
close(gcf);

function Blank_Output_MS_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function Blank_Output_MS_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Data_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function Data_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
