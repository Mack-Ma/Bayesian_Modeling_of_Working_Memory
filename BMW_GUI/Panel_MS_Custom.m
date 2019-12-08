%% Panel_MC_Custom (Model Specification: Custom Models)
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

function varargout = Panel_MS_Custom(varargin)

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Panel_MS_Custom_OpeningFcn, ...
                   'gui_OutputFcn',  @Panel_MS_Custom_OutputFcn, ...
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


% --- Executes just before Panel_MS_Custom is made visible.
function Panel_MS_Custom_OpeningFcn(hObject, ~, handles, varargin)

handles.output = hObject;
guidata(hObject, handles);

WindowSize=get(gcf,'Position');
WindowSize(1:2)=[39.8 10];
set(gcf,'Position',WindowSize)

global bmw AllModelName;
cd([bmw.path, 'BMW_Models\Custom_Models'])
AllModel=dir('*.m');
AllModelName={AllModel.name};
set(handles.List_Model_Custom,'String',AllModelName)


% --- Outputs from this function are returned to the command line.
function varargout = Panel_MS_Custom_OutputFcn(~, ~, handles) 

varargout{1} = handles.output;

% --- Executes on selection change in List_Model_Custom.
function List_Model_Custom_Callback(~, ~, ~)

% --- Executes during object creation, after setting all properties.
function List_Model_Custom_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(~, ~, handles)
global AllModelName MOI
MOI.Model=AllModelName{get(handles.List_Model_Custom,'Value')};
close(gcf);
