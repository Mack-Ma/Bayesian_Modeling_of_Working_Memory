function varargout = BMW_Main(varargin)

% Edit the above text to modify the response to help BMW_Main

% Last Modified by GUIDE v2.5 28-Mar-2020 21:21:13

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @BMW_Main_OpeningFcn, ...
                   'gui_OutputFcn',  @BMW_Main_OutputFcn, ...
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


% --- Executes just before BMW_Main is made visible.
function BMW_Main_OpeningFcn(hObject, eventdata, handles, varargin)

% Choose default command line output for BMW_Main
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = BMW_Main_OutputFcn(hObject, eventdata, handles) 
% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in Fit_Compare_BMW.
function Fit_Compare_BMW_Callback(hObject, eventdata, handles)
BMW_Run;

% --- Executes on button press in Change_Settings_BMW.
function Change_Settings_BMW_Callback(hObject, eventdata, handles)
BMW_Change_Settings;

% --- Executes on button press in About_BMW.
function About_BMW_Callback(hObject, eventdata, handles)
BMW_About;

% --- Executes on button press in Manual_BMW.
function Manual_BMW_Callback(hObject, eventdata, handles)
BMW('Manual');
