%% Panel_MC_CR (Model Specification: Continuous Recall)
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

function varargout = Panel_MS_CR(varargin)

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Panel_MS_CR_OpeningFcn, ...
                   'gui_OutputFcn',  @Panel_MS_CR_OutputFcn, ...
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

% --- Executes just before Panel_MS_CR is made visible.
function Panel_MS_CR_OpeningFcn(hObject, ~, handles, varargin)

global AllModelName
AllModelName={'Item Limit','Standard Mixture','Equal Precision','Slots-plus-Averaging',...
    'Variable Precision','Variable Precision with Capacity', 'Categorical Slots-plus-Averaging',...
    'Categorical Variable Precision', 'Categorical Variable Precision with Capacity'};

handles.output = hObject;
guidata(hObject, handles);

WindowSize=get(gcf,'Position');
WindowSize(1:2)=[39.8 10];
set(gcf,'Position',WindowSize)

set(handles.List_Model_CR,'String',AllModelName)

% --- Outputs from this function are returned to the command line.
function varargout = Panel_MS_CR_OutputFcn(~, ~, handles) 

varargout{1} = handles.output;

% --- Executes on selection change in List_Model_CR.
function List_Model_CR_Callback(~, ~, ~)

% --- Executes during object creation, after setting all properties.
function List_Model_CR_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(~, ~, handles)
global AllModelName Derivatives MOI
MOI.Derivatives=Derivatives;
MOI.Model=AllModelName{get(handles.List_Model_CR,'Value')};
close(gcf);

% --- Executes on button press in Deriv_Swap.
function Deriv_Swap_Callback(hObject, ~, ~)
global Derivatives
if hObject.Value==1
    Derivatives.Swap=1;
end

% --- Executes on button press in Deriv_Bias.
function Deriv_Bias_Callback(hObject, ~, ~)
global Derivatives
if hObject.Value==1
    Derivatives.Bias=1;
end

% --- Executes on button press in Deriv_BiasF.
function Deriv_BiasF_Callback(hObject, ~, ~)
global Derivatives
if hObject.Value==1
    Derivatives.BiasF=1;
end

% --- Executes on button press in Deriv_PrecF.
function Deriv_PrecF_Callback(hObject, ~, ~)
global Derivatives
if hObject.Value==1
    Derivatives.PrecF=1;
end
