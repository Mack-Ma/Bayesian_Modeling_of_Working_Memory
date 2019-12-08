%% Panel_Conf (Configuration: Main Panel)
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

function varargout = Panel_Conf(varargin)

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Panel_Conf_OpeningFcn, ...
                   'gui_OutputFcn',  @Panel_Conf_OutputFcn, ...
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


% --- Executes just before Panel_Conf is made visible.
function Panel_Conf_OpeningFcn(hObject, eventdata, handles, varargin)

handles.output = hObject;
guidata(hObject, handles);

WindowSize=get(gcf,'Position');
WindowSize(1:2)=[150 25];
set(gcf,'Position',WindowSize)

AllCriteria={'Default','Define'};
AllConstraints={'Default','Define'};
AllFitOptions={'Default','Define'};
set(handles.List_Criteria,'String',AllCriteria)
set(handles.List_Constraints,'String',AllConstraints)
set(handles.List_FitOptions,'String',AllFitOptions)

% --- Outputs from this function are returned to the command line.
function varargout = Panel_Conf_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

% --- Executes on selection change in List_Criteria.
function List_Criteria_Callback(hObject, eventdata, handles)
global Criteria
switch hObject.Value
    case 1
        Criteria.default=1;
    case 2
        Criteria.default=0;
        Panel_Conf_Criteria;
end

% --- Executes during object creation, after setting all properties.
function List_Criteria_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Done_MS.
function Done_MS_Callback(hObject, eventdata, handles)
% hObject    handle to Done_MS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Criteria Constraints FitOptions
Conf_BMW.InputFile=get(handles.Blank_File_Conf,'String');
Conf_BMW.Criteria=Criteria;
Conf_BMW.Constraints=Constraints;
Conf_BMW.FitOptions=FitOptions;
Configuration_BMW(Conf_BMW);
close(gcf);

function Blank_File_Conf_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function Blank_File_Conf_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in List_FitOptions.
function List_FitOptions_Callback(hObject, eventdata, handles)
global FitOptions
switch hObject.Value
    case 1
        FitOptions.default=1;
    case 2
        FitOptions.default=0;
        Panel_Conf_FitOptions;
end

% --- Executes during object creation, after setting all properties.
function List_FitOptions_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in List_Constraints.
function List_Constraints_Callback(hObject, eventdata, handles)
global Constraints
switch hObject.Value
    case 1
        Constraints.default=1;
    case 2
        Constraints.default=0;
        Panel_Conf_Constraints;
end

% --- Executes during object creation, after setting all properties.
function List_Constraints_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
