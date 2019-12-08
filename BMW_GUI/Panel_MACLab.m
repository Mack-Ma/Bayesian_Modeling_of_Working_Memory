%% Panel_MACLab (Lab Intro)
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

function varargout = Panel_MACLab(varargin)

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Panel_MACLab_OpeningFcn, ...
                   'gui_OutputFcn',  @Panel_MACLab_OutputFcn, ...
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


% --- Executes just before Panel_MACLab is made visible.
function Panel_MACLab_OpeningFcn(hObject, eventdata, handles, varargin)

handles.output = hObject;
guidata(hObject, handles);

WindowSize=get(gcf,'Position');
WindowSize(1:2)=[150 25];
set(gcf,'Position',WindowSize)

% --- Outputs from this function are returned to the command line.
function varargout = Panel_MACLab_OutputFcn(hObject, eventdata, handles) 

varargout{1} = handles.output;
