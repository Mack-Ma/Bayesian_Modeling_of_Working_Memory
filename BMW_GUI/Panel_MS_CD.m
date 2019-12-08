%% Panel_MC_CD (Model Specification: Change Detection)
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

function varargout = Panel_MS_CD(varargin)

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Panel_MS_CD_OpeningFcn, ...
                   'gui_OutputFcn',  @Panel_MS_CD_OutputFcn, ...
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


% --- Executes just before Panel_MS_CD is made visible.
function Panel_MS_CD_OpeningFcn(hObject, ~, handles, varargin)

handles.output = hObject;
guidata(hObject, handles);

WindowSize=get(gcf,'Position');
WindowSize(1:2)=[39.8 10];
set(gcf,'Position',WindowSize)

global AllModel
AllModel={'Fixed-Capacity (Full-Display)', 'Fixed-Capacity (Single-Probe)', 'Fixed-Capacity (Central Probe)', 'Signal Detection'};
set(handles.List_Model_CD,'String',AllModel)

% --- Outputs from this function are returned to the command line.
function varargout = Panel_MS_CD_OutputFcn(~, ~, handles) 

varargout{1} = handles.output;

% --- Executes on selection change in List_Model_CD.
function List_Model_CD_Callback(~, ~, ~)

% --- Executes during object creation, after setting all properties.
function List_Model_CD_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(~, ~, handles)
global Derivatives AllModel MOI
MOI.Derivatives=Derivatives;
MOI.Model=AllModel{get(handles.List_Model_CD,'Value')};
close(gcf);

% --- Executes on button press in Deriv_Ensem.
function Deriv_Ensem_Callback(hObject, ~, ~)
global Derivatives
if hObject.Value==1
    Derivatives.Ensemble=1;
end

% --- Executes on button press in Deriv_Lapse.
function Deriv_Lapse_Callback(hObject, ~, ~)
global Derivatives
if hObject.Value==1
    Derivatives.Lapse=1;
end
