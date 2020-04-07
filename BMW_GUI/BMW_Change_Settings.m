function varargout = BMW_Change_Settings(varargin)

% Last Modified by GUIDE v2.5 29-Mar-2020 02:07:11

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @BMW_Change_Settings_OpeningFcn, ...
                   'gui_OutputFcn',  @BMW_Change_Settings_OutputFcn, ...
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


% --- Executes just before BMW_Change_Settings is made visible.
function BMW_Change_Settings_OpeningFcn(hObject, eventdata, handles, varargin)

% Choose default command line output for BMW_Change_Settings
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes BMW_Change_Settings wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = BMW_Change_Settings_OutputFcn(hObject, eventdata, handles) 

% Get default command line output from handles structure
varargout{1} = handles.output;


function MCMCO_NBB_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function MCMCO_NBB_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function MCMCO_NBS_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function MCMCO_NBS_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function MCMCO_Chain_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function MCMCO_Chain_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function MCMCO_NS_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function MCMCO_NS_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function MCMCO_R_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function MCMCO_R_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in popupmenu5.
function popupmenu5_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function popupmenu5_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in RFXBMS.
function RFXBMS_Callback(hObject, eventdata, handles)

% --- Executes on button press in checkbox58.
function checkbox58_Callback(hObject, eventdata, handles)

% --- Executes on button press in checkbox65.
function checkbox65_Callback(hObject, eventdata, handles)

% --- Executes on button press in checkbox64.
function checkbox64_Callback(hObject, eventdata, handles)

% --- Executes on button press in checkbox63.
function checkbox63_Callback(hObject, eventdata, handles)

% --- Executes on button press in checkbox62.
function checkbox62_Callback(hObject, eventdata, handles)

% --- Executes on button press in checkbox61.
function checkbox61_Callback(hObject, eventdata, handles)

% --- Executes on button press in checkbox60.
function checkbox60_Callback(hObject, eventdata, handles)

% --- Executes on button press in checkbox59.
function checkbox59_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox59 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox59


% --- Executes on button press in checkbox57.
function checkbox57_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox57 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox57


% --- Executes on selection change in FO_Algorithm.
function FO_Algorithm_Callback(hObject, eventdata, handles)
% hObject    handle to FO_Algorithm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns FO_Algorithm contents as cell array
%        contents{get(hObject,'Value')} returns selected item from FO_Algorithm


% --- Executes during object creation, after setting all properties.
function FO_Algorithm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FO_Algorithm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu3.
function popupmenu3_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu3


% --- Executes during object creation, after setting all properties.
function popupmenu3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu4.
function popupmenu4_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu4 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu4


% --- Executes during object creation, after setting all properties.
function popupmenu4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in MV_CD_Lapse.
function MV_CD_Lapse_Callback(hObject, eventdata, handles)
% hObject    handle to MV_CD_Lapse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of MV_CD_Lapse


% --- Executes on button press in MV_CD_Ensemble.
function MV_CD_Ensemble_Callback(hObject, eventdata, handles)
% hObject    handle to MV_CD_Ensemble (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of MV_CD_Ensemble


% --- Executes on button press in MV_CD_None.
function MV_CD_None_Callback(hObject, eventdata, handles)
% hObject    handle to MV_CD_None (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of MV_CD_None


% --- Executes on button press in Fixed_Capacity_Full_Display.
function Fixed_Capacity_Full_Display_Callback(hObject, eventdata, handles)
% hObject    handle to Fixed_Capacity_Full_Display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Fixed_Capacity_Full_Display


% --- Executes on button press in Fixed_Capacity_Single_Probe.
function Fixed_Capacity_Single_Probe_Callback(hObject, eventdata, handles)
% hObject    handle to Fixed_Capacity_Single_Probe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Fixed_Capacity_Single_Probe


% --- Executes on button press in checkbox28.
function checkbox28_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox28 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox28


% --- Executes on button press in Signal_Detection.
function Signal_Detection_Callback(hObject, eventdata, handles)
% hObject    handle to Signal_Detection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Signal_Detection


% --- Executes on button press in MV_CR_Bias.
function MV_CR_Bias_Callback(hObject, eventdata, handles)

% --- Executes on button press in MV_CR_BiasF.
function MV_CR_BiasF_Callback(hObject, eventdata, handles)

% --- Executes on button press in MV_CR_PrecF.
function MV_CR_PrecF_Callback(hObject, eventdata, handles)

% --- Executes on button press in MV_CR_None.
function MV_CR_None_Callback(hObject, eventdata, handles)

% --- Executes on button press in MV_CR_Swap.
function MV_CR_Swap_Callback(hObject, eventdata, handles)

% --- Executes on button press in Item_Limit.
function Item_Limit_Callback(hObject, eventdata, handles)

% --- Executes on button press in Standard_Mixture.
function Standard_Mixture_Callback(hObject, eventdata, handles)

% --- Executes on button press in Slots_plus_Averaging.
function Slots_plus_Averaging_Callback(hObject, eventdata, handles)

% --- Executes on button press in checkbox4.
function checkbox4_Callback(hObject, eventdata, handles)

% --- Executes on button press in Variable_Precision.
function Variable_Precision_Callback(hObject, eventdata, handles)

% --- Executes on button press in Variable_Precision_with_Capacity.
function Variable_Precision_with_Capacity_Callback(hObject, eventdata, handles)

% --- Executes on button press in Categorical_Slots_plus_Averaging_BV.
function Categorical_Slots_plus_Averaging_BV_Callback(hObject, eventdata, handles)

% --- Executes on button press in Categorical_Variable_Precision.
function Categorical_Variable_Precision_Callback(hObject, eventdata, handles)

% --- Executes on button press in Categorical_Variable_Precision_with_Capacity.
function Categorical_Variable_Precision_with_Capacity_Callback(hObject, eventdata, handles)
