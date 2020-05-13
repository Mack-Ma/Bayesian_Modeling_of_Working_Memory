function varargout = BMW_Run(varargin)

% Last Modified by GUIDE v2.5 27-Apr-2020 01:48:48

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @BMW_Run_OpeningFcn, ...
                   'gui_OutputFcn',  @BMW_Run_OutputFcn, ...
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


% --- Executes just before BMW_Run is made visible.
function BMW_Run_OpeningFcn(hObject, eventdata, handles, varargin)

% Choose default command line output for BMW_Run
handles.output = hObject;

global BMW_Preferences
% Load settings
load('BMW_Preferences')

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes BMW_Run wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = BMW_Run_OutputFcn(hObject, eventdata, handles) 

% Get default command line output from handles structure
varargout{1} = handles.output;

function InputDir_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function InputDir_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function OutputDir_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function OutputDir_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in Fit.
function Fit_Callback(hObject, eventdata, handles)
global BMW_Preferences
BMW_Input=BMW_Preferences;
BMW_Input.Data=get(handles.InputDir,'String');
MA_BMW=cell(length(BMW_Input.ModelSpace),BMW_Input.FitOptions.Repeat);
for j=1:length(BMW_Input.ModelSpace)
    for i=1:BMW_Input.FitOptions.Repeat
        fprintf('Model: %s; Repeat: %d',BMW_Input.ModelSpace{j},i)
        BMW_Input.Model.Model=BMW_Input.ModelSpace{j};
        MA_BMW{i,j}=ModelFit_BMW(Configuration_BMW(BMW_Input));
    end
end
CurDir=cd;
cd(get(handles.OutputDir,'String'));
save('MA_BMW','MA_BMW');
cd(CurDir);

% --- Executes on button press in Compare.
function Compare_Callback(hObject, eventdata, handles)
BMW_Input=get(handles.InputDir,'String');
MC_BMW=ModelComparison_BMW(BMW_Input);
CurDir=cd;
cd(get(handles.OutputDir,'String'));
save('MC_BMW','MC_BMW');
cd(CurDir);
