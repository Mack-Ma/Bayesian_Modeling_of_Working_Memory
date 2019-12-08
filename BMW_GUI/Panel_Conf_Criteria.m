%% Panel_Conf_Criteria (Configuration: Criteria)
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

function varargout = Panel_Conf_Criteria(varargin)

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Panel_Conf_Criteria_OpeningFcn, ...
                   'gui_OutputFcn',  @Panel_Conf_Criteria_OutputFcn, ...
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

% --- Executes just before Panel_Conf_Criteria is made visible.
function Panel_Conf_Criteria_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

WindowSize=get(gcf,'Position');
WindowSize(1:2)=[39.8 20];
set(gcf,'Position',WindowSize)

% --- Outputs from this function are returned to the command line.
function varargout = Panel_Conf_Criteria_OutputFcn(hObject, eventdata, handles) 

varargout{1} = handles.output;

% --- Executes on button press in LLH.
function LLH_Callback(hObject, eventdata, handles)

% --- Executes on button press in AIC.
function AIC_Callback(hObject, eventdata, handles)

% --- Executes on button press in AICc.
function AICc_Callback(hObject, eventdata, handles)

% --- Executes on button press in BIC.
function BIC_Callback(hObject, eventdata, handles)

% --- Executes on button press in DIC.
function DIC_Callback(hObject, eventdata, handles)

% --- Executes on button press in Done_Conf_Criteria.
function Done_Conf_Criteria_Callback(hObject, eventdata, handles)
global Criteria
Criteria.LME=get(handles.LME,'Value');
Criteria.LLH=get(handles.LLH,'Value');
Criteria.AIC=get(handles.AIC,'Value');
Criteria.AICc=get(handles.AICc,'Value');
Criteria.BIC=get(handles.BIC,'Value');
Criteria.DIC=get(handles.DIC,'Value');
Criteria.WAIC=get(handles.WAIC,'Value');
close(gcf);

% --- Executes on button press in WAIC.
function WAIC_Callback(hObject, eventdata, handles)

% --- Executes on button press in LME.
function LME_Callback(hObject, eventdata, handles)
