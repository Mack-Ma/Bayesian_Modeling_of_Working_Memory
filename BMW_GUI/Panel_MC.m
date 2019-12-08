%% Panel_MC (Model Comparison)
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

function varargout = Panel_MC(varargin)

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Panel_MC_OpeningFcn, ...
                   'gui_OutputFcn',  @Panel_MC_OutputFcn, ...
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


% --- Executes just before Panel_MC is made visible.
function Panel_MC_OpeningFcn(hObject, eventdata, handles, varargin)

handles.output = hObject;
guidata(hObject, handles);

WindowSize=get(gcf,'Position');
WindowSize(1:2)=[150 25];
set(gcf,'Position',WindowSize)

global AllCriteria
AllModelNames={'Default','Define'};
AllCriteria={'LLH','AIC','AICc','BIC','DIC','WAIC','LME'};
set(handles.Model_Names,'String',AllModelNames)
set(handles.List_Criterion,'String',AllCriteria)
% UIWAIT makes Panel_MC wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Panel_MC_OutputFcn(hObject, eventdata, handles) 
% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in Done_MC.
function Done_MC_Callback(hObject, eventdata, handles)
global Model_Space Model_Names AllCriteria
MC_BMW.MS=Model_Space;
MC_BMW.MN=Model_Names;
MC_BMW.Criterion=AllCriteria{get(handles.List_Criterion,'Value')};
MC_BMW.OutputDir=get(handles.OutputDir,'String');
ModelComparison_BMW(MC_BMW);
close(gcf);

function Model_Space_Callback(hObject, eventdata, handles)

Panel_MC_MS;

% --- Executes during object creation, after setting all properties.
function Model_Space_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function List_Criterion_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function List_Criterion_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function OutputDir_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function OutputDir_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function Model_Names_Callback(hObject, eventdata, handles)

global Model_Names
if strcmp(hObject.String,'Define')
    Panel_MC_ModelNames;
elseif strcmp(hObject.String,'Default')
    Model_Names=0;
end

% --- Executes during object creation, after setting all properties.
function Model_Names_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
