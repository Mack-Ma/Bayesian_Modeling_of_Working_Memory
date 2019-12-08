%% Panel_Main (Main Panel)
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

function varargout = Panel_Main(varargin)

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Panel_Main_OpeningFcn, ...
                   'gui_OutputFcn',  @Panel_Main_OutputFcn, ...
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

% --- Executes just before Panel_Main is made visible.
function Panel_Main_OpeningFcn(hObject, eventdata, handles, varargin)

handles.output = hObject;
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = Panel_Main_OutputFcn(hObject, eventdata, handles) 

varargout{1} = handles.output;

% --- Executes on button press in MD.
function MD_Callback(hObject, eventdata, handles)
Panel_MS;

% --- Executes on button press in Conf.
function Conf_Callback(hObject, eventdata, handles)
Panel_Conf;

% --- Executes on button press in MA.
function MA_Callback(hObject, eventdata, handles)
Panel_MA;

% --- Executes on button press in MC.
function MC_Callback(hObject, eventdata, handles)
Panel_MC;

% --- Executes on button press in Vis.
function Vis_Callback(hObject, eventdata, handles)
Panel_Vis;

% --- Executes on button press in MAC_Lab.
function MAC_Lab_Callback(hObject, eventdata, handles)
Panel_MACLab;
