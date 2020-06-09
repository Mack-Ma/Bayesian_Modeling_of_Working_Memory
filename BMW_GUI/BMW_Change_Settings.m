function varargout = BMW_Change_Settings(varargin)

% Last Modified by GUIDE v2.5 08-Jun-2020 00:15:54

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

global BMW_Preferences
load('BMW_Preferences.mat');
BMWSet(BMW_Preferences,handles)

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes BMW_Change_Settings wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = BMW_Change_Settings_OutputFcn(hObject, eventdata, handles)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in RFXBMS.
function RFXBMS_Callback(hObject, eventdata, handles)

% --- Executes on button press in AIC.
function AIC_Callback(hObject, eventdata, handles)

% --- Executes on button press in LME_GHM.
function LME_GHM_Callback(hObject, eventdata, handles)

% --- Executes on button press in WAIC2.
function WAIC2_Callback(hObject, eventdata, handles)

% --- Executes on button press in WAIC1.
function WAIC1_Callback(hObject, eventdata, handles)

% --- Executes on button press in DICs.
function DICs_Callback(hObject, eventdata, handles)

% --- Executes on button press in BIC.
function BIC_Callback(hObject, eventdata, handles)

% --- Executes on button press in DIC1.
function DIC1_Callback(hObject, eventdata, handles)

% --- Executes on button press in AICc.
function AICc_Callback(hObject, eventdata, handles)

% --- Executes on button press in LLH.
function LLH_Callback(hObject, eventdata, handles)

% --- Executes on selection change in FitOptions_Algorithm.
function FitOptions_Algorithm_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function FitOptions_Algorithm_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in FitOptions_Method.
function FitOptions_Method_Callback(hObject, eventdata, handles)
switch hObject.Value
    case 1
        set(handles.FitOptions_Prior,'Enable','on')
    case 2
        set(handles.FitOptions_Prior,'Enable','off')
end

% --- Executes during object creation, after setting all properties.
function FitOptions_Method_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in FitOptions_Verbosity.
function FitOptions_Verbosity_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function FitOptions_Verbosity_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function FitOptions_Repeat_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function FitOptions_Repeat_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in Lapse.
function Lapse_Callback(hObject, eventdata, handles)

% --- Executes on button press in Ensemble.
function Ensemble_Callback(hObject, eventdata, handles)

% --- Executes on button press in Fixed_Capacity_Full_Display.
function Fixed_Capacity_Full_Display_Callback(hObject, eventdata, handles)

% --- Executes on button press in Fixed_Capacity_Single_Probe.
function Fixed_Capacity_Single_Probe_Callback(hObject, eventdata, handles)

% --- Executes on button press in Signal_Detection.
function Signal_Detection_Callback(hObject, eventdata, handles)

% --- Executes on button press in Bias.
function Bias_Callback(hObject, eventdata, handles)

% --- Executes on button press in BiasF.
function BiasF_Callback(hObject, eventdata, handles)

% --- Executes on button press in PrecF.
function PrecF_Callback(hObject, eventdata, handles)

% --- Executes on button press in Swap.
function Swap_Callback(hObject, eventdata, handles)

% --- Executes on button press in Item_Limit.
function Item_Limit_Callback(hObject, eventdata, handles)

% --- Executes on button press in Standard_Mixture.
function Standard_Mixture_Callback(hObject, eventdata, handles)

% --- Executes on button press in Slots_plus_Averaging.
function Slots_plus_Averaging_Callback(hObject, eventdata, handles)

% --- Executes on button press in Equal_Precision.
function Equal_Precision_Callback(hObject, eventdata, handles)

% --- Executes on button press in Variable_Precision.
function Variable_Precision_Callback(hObject, eventdata, handles)

% --- Executes on button press in Variable_Precision_with_Capacity.
function Variable_Precision_with_Capacity_Callback(hObject, eventdata, handles)
% --- Executes on button press in ContinuousK.

% --- Executes on button press in Category_Only_with_Capacity.
function Category_Only_with_Capacity_Callback(hObject, eventdata, handles)

function ContinuousK_Callback(hObject, eventdata, handles)

% --- Executes on button press in ResponseNoise.
function ResponseNoise_Callback(hObject, eventdata, handles)

% --- Executes on button press in Category_WI.
function Category_WI_Callback(hObject, eventdata, handles)

% --- Executes on button press in Category_BI.
function Category_BI_Callback(hObject, eventdata, handles)

% --- Executes on button press in Category_Only.
function Category_Only_Callback(hObject, eventdata, handles)

% --- Executes on button press in Fixed_Capacity_Central_Probe.
function Fixed_Capacity_Central_Probe_Callback(hObject, eventdata, handles)

% --- Executes on button press in DIC2.
function DIC2_Callback(hObject, eventdata, handles)

% --- Executes on button press in LME_BS.
function LME_BS_Callback(hObject, eventdata, handles)

% --- Executes on selection change in MCMC_GR.
function MCMC_GR_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function MCMC_GR_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function MCMC_Sample_Manual_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function MCMC_Sample_Manual_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in MCMC_Tol_Mode.
function MCMC_Tol_Mode_Callback(hObject, eventdata, handles)
switch hObject.Value
    case 1
        set(handles.MCMC_Tol_Manual,'String','');
        set(handles.MCMC_Tol_Manual,'Enable','off');
    case 2
        set(handles.MCMC_Tol_Manual,'Enable','on');
end

% --- Executes during object creation, after setting all properties.
function MCMC_Tol_Mode_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in MCMC_Sample_Mode.
function MCMC_Sample_Mode_Callback(hObject, eventdata, handles)
switch hObject.Value
    case 1
        set(handles.MCMC_Sample_Manual,'String','');
        set(handles.MCMC_Sample_Manual,'Enable','off');
    case 2
        set(handles.MCMC_Sample_Manual,'Enable','on');
end

% --- Executes during object creation, after setting all properties.
function MCMC_Sample_Mode_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function MCMC_Tol_Manual_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function MCMC_Tol_Manual_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function MCMC_Ncore_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function MCMC_Ncore_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in MCMC_TransformMethod.
function MCMC_TransformMethod_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function MCMC_TransformMethod_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function MCMC_BurninBatches_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function MCMC_BurninBatches_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function MCMC_BurninSample_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function MCMC_BurninSample_Callback(hObject, eventdata, handles)

% --- Executes on selection change in FitOptions_Prior.
function FitOptions_Prior_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function FitOptions_Prior_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in Change.
function Change_Callback(hObject, eventdata, handles)
BMW_Preferences=BMWGet(handles);
CurDir=cd;
CSDir=which('BMW_Preferences.mat');
CSDir=CSDir(1:end-length('BMW_Preferences.mat'));
cd(CSDir);
save('BMW_Preferences','BMW_Preferences');
cd(CurDir);
close

% --- Executes on button press in Defaults.
function Defaults_Callback(hObject, eventdata, handles)
BMW_Preferences.ModelSpace={'Standard Mixture'};
BMW_Preferences.Model.Variants={};
BMW_Preferences.FitOptions.Method='MAP';
BMW_Preferences.FitOptions.UniformPrior=0;
BMW_Preferences.FitOptions.Repeat=1;
BMW_Preferences.FitOptions.MCMCOptions.Transform='Probit';
BMW_Preferences.FitOptions.Display='iter';
BMW_Preferences.FitOptions.MCMCOptions.Convergence.Diagnostic='GR'; 
BMW_Preferences.FitOptions.MCMCOptions.Convergence.Nbatchburnin=200; 
BMW_Preferences.FitOptions.MCMCOptions.Convergence.Nmaxbatchburnin=50; 
BMW_Preferences.FitOptions.MCMCOptions.Ncore=4;
BMW_Preferences.FitOptions.Algorithm='DE-MCMC';
BMW_Preferences.Criteria={'DIC2','WAIC2','LME_GHM'};
BMW_Preferences.FitOptions.RFXBMS=1;
BMWSet(BMW_Preferences,handles);
CurDir=cd;
CSDir=which('BMW_Preferences.mat');
CSDir=CSDir(1:end-length('BMW_Preferences.mat'));
cd(CSDir);
save('BMW_Preferences','BMW_Preferences');
cd(CurDir);

function BMWSet(BMW_Preferences,handles)
% checkboxes
set(handles.RFXBMS,'Value',BMW_Preferences.FitOptions.RFXBMS)
set(handles.LLH,'Value',any(strcmp(BMW_Preferences.Criteria,'LLH')))
set(handles.DIC1,'Value',any(strcmp(BMW_Preferences.Criteria,'DIC1')))
set(handles.DIC2,'Value',any(strcmp(BMW_Preferences.Criteria,'DIC2')))
set(handles.DICs,'Value',any(strcmp(BMW_Preferences.Criteria,'DICs')))
set(handles.AIC,'Value',any(strcmp(BMW_Preferences.Criteria,'AIC')))
set(handles.AICc,'Value',any(strcmp(BMW_Preferences.Criteria,'AICc')))
set(handles.BIC,'Value',any(strcmp(BMW_Preferences.Criteria,'BIC')))
set(handles.WAIC1,'Value',any(strcmp(BMW_Preferences.Criteria,'WAIC1')))
set(handles.WAIC2,'Value',any(strcmp(BMW_Preferences.Criteria,'WAIC2')))
set(handles.LME_BS,'Value',any(strcmp(BMW_Preferences.Criteria,'LME_BS')))
set(handles.LME_GHM,'Value',any(strcmp(BMW_Preferences.Criteria,'LME_GHM')))
set(handles.Item_Limit,'Value',any(strcmp(BMW_Preferences.ModelSpace,'Item Limit')))
set(handles.Standard_Mixture,'Value',any(strcmp(BMW_Preferences.ModelSpace,'Standard Mixture')))
set(handles.Slots_plus_Averaging,'Value',any(strcmp(BMW_Preferences.ModelSpace,'Slots-plus-Averaging')))
set(handles.Equal_Precision,'Value',any(strcmp(BMW_Preferences.ModelSpace,'Equal Precision')))
set(handles.Variable_Precision,'Value',any(strcmp(BMW_Preferences.ModelSpace,'Variable Precision')))
set(handles.Category_Only,'Value',any(strcmp(BMW_Preferences.ModelSpace,'Category-Only')))
set(handles.Category_Only_with_Capacity,'Value',any(strcmp(BMW_Preferences.ModelSpace,'Category-Only (with Capacity)')))
set(handles.Variable_Precision_with_Capacity,'Value',any(strcmp(BMW_Preferences.ModelSpace,'Variable Precision with Capacity')))
set(handles.Fixed_Capacity_Central_Probe,'Value',any(strcmp(BMW_Preferences.ModelSpace,'Fixed Capacity (Central-Probe)')))
set(handles.Fixed_Capacity_Single_Probe,'Value',any(strcmp(BMW_Preferences.ModelSpace,'Fixed Capacity (Single-Probe)')))
set(handles.Fixed_Capacity_Full_Display,'Value',any(strcmp(BMW_Preferences.ModelSpace,'Fixed Capacity (Full-Display)')))
set(handles.Signal_Detection,'Value',any(strcmp(BMW_Preferences.ModelSpace,'Signal Detection')))
set(handles.Bias,'Value',any(strcmp(BMW_Preferences.Model.Variants,'Bias')))
set(handles.BiasF,'Value',any(strcmp(BMW_Preferences.Model.Variants,'BiasF')))
set(handles.PrecF,'Value',any(strcmp(BMW_Preferences.Model.Variants,'PrecF')))
set(handles.Swap,'Value',any(strcmp(BMW_Preferences.Model.Variants,'Swap')))
set(handles.ContinuousK,'Value',any(strcmp(BMW_Preferences.Model.Variants,'ContinuousK')))
set(handles.ResponseNoise,'Value',any(strcmp(BMW_Preferences.Model.Variants,'ResponseNoise')))
set(handles.Category_WI,'Value',any(strcmp(BMW_Preferences.Model.Variants,'Category (Within-Item)')))
set(handles.Category_BI,'Value',any(strcmp(BMW_Preferences.Model.Variants,'Category (Between-Item)')))
set(handles.Lapse,'Value',any(strcmp(BMW_Preferences.Model.Variants,'Lapse')))
set(handles.Ensemble,'Value',any(strcmp(BMW_Preferences.Model.Variants,'Ensemble')))
% mcmc options
switch BMW_Preferences.FitOptions.MCMCOptions.Convergence.Diagnostic
    case 'GR'
        ind_GR=1;
    case 'GRL'
        ind_GR=2;
end
set(handles.MCMC_GR,'Value',ind_GR)
if isfield(BMW_Preferences.FitOptions.MCMCOptions,'Nsample')
    MCMC_Nsample=BMW_Preferences.FitOptions.MCMCOptions.Nsample;
    set(handles.MCMC_Sample_Manual,'Enable','on')
    ind_Nsample=2;
else
    MCMC_Nsample='';
    set(handles.MCMC_Sample_Manual,'Enable','off')
    ind_Nsample=1;
end
set(handles.MCMC_Sample_Mode,'Value',ind_Nsample)
set(handles.MCMC_Sample_Manual,'String',num2str(MCMC_Nsample))
if isfield(BMW_Preferences.FitOptions.MCMCOptions.Convergence,'Tol')
    MCMC_Tol=BMW_Preferences.FitOptions.MCMCOptions.Tol;
    set(handles.MCMC_Tol_Manual,'Enable','on')
    ind_Tol=2;
else
    MCMC_Tol='';
    set(handles.MCMC_Tol_Manual,'Enable','off')
    ind_Tol=1;
end
set(handles.MCMC_Tol_Mode,'Value',ind_Tol)
set(handles.MCMC_Tol_Manual,'String',num2str(MCMC_Tol))
switch BMW_Preferences.FitOptions.MCMCOptions.Transform
    case 'Probit'
        MCMC_Transform=1;
    case 'Logit'
        MCMC_Transform=2;
    case 'Fisher'
        MCMC_Transform=3;
    case 'NoTransform'
        MCMC_Transform=4;
end
set(handles.MCMC_TransformMethod,'Value',MCMC_Transform)
set(handles.MCMC_Ncore,'String',num2str(BMW_Preferences.FitOptions.MCMCOptions.Ncore))
set(handles.MCMC_BurninBatches,'String',num2str(BMW_Preferences.FitOptions.MCMCOptions.Convergence.Nbatchburnin))
set(handles.MCMC_BurninSample,'String',num2str(BMW_Preferences.FitOptions.MCMCOptions.Convergence.Nmaxbatchburnin))
% fit options
set(handles.FitOptions_Repeat,'String',num2str(BMW_Preferences.FitOptions.Repeat))
set(handles.FitOptions_Prior,'Value',BMW_Preferences.FitOptions.UniformPrior+1)
switch BMW_Preferences.FitOptions.Method
    case 'MAP'
        ind_Method=1;
        set(handles.FitOptions_Prior,'Enable','on')
    case 'MLE'
        ind_Method=2;
        set(handles.FitOptions_Prior,'Enable','off')
end
set(handles.FitOptions_Method,'Value',ind_Method)
switch BMW_Preferences.FitOptions.Display
    case 'iter'
        ind_Verbosity=1;
    case 'off'
        ind_Verbosity=2;
    case 'detail'
        ind_Verbosity=3;
end
set(handles.FitOptions_Verbosity,'Value',ind_Verbosity)

function BMW_Preferences=BMWGet(handles)
% mcmc options
switch get(handles.MCMC_GR,'Value')
    case 1
        BMW_Preferences.FitOptions.MCMCOptions.Convergence.Diagnostic='GR';
    case 2
        BMW_Preferences.FitOptions.MCMCOptions.Convergence.Diagnostic='GRL';
end
switch get(handles.MCMC_Sample_Mode,'Value')
    case 2
        BMW_Preferences.FitOptions.MCMCOptions.Convergence.Nsample=str2double(get(handles.MCMC_Sample_Manual,'String'));
end
switch get(handles.MCMC_Tol_Mode,'Value')
    case 2
        BMW_Preferences.FitOptions.MCMCOptions.Convergence.Tol=str2double(get(handles.MCMC_Tol_Manual,'String'));
end
switch get(handles.MCMC_TransformMethod,'Value')
    case 1
        BMW_Preferences.FitOptions.MCMCOptions.Transform='Probit';
    case 2
        BMW_Preferences.FitOptions.MCMCOptions.Transform='Logit';
    case 3
        BMW_Preferences.FitOptions.MCMCOptions.Transform='Fisher';
    case 4
        BMW_Preferences.FitOptions.MCMCOptions.Transform='NoTransform';
end
BMW_Preferences.FitOptions.MCMCOptions.Convergence.Nbatchburnin=str2double(get(handles.MCMC_BurninSample,'String'));
BMW_Preferences.FitOptions.MCMCOptions.Convergence.Nmaxbatchburnin=str2double(get(handles.MCMC_BurninBatches,'String'));
BMW_Preferences.FitOptions.MCMCOptions.Ncore=str2double(get(handles.MCMC_Ncore,'String'));
% fit options
BMW_Preferences.FitOptions.RFXBMS=get(handles.RFXBMS,'Value');
BMW_Preferences.FitOptions.Repeat=str2double(get(handles.FitOptions_Repeat,'String'));
switch get(handles.FitOptions_Method,'Value')
    case 1
        BMW_Preferences.FitOptions.Method='MAP';
    case 2
        BMW_Preferences.FitOptions.Method='MLE';
end
BMW_Preferences.FitOptions.UniformPrior=get(handles.FitOptions_Prior,'Value')-1;
switch get(handles.FitOptions_Verbosity,'Value')
    case 1
        BMW_Preferences.FitOptions.Display='iter';
    case 2
        BMW_Preferences.FitOptions.Display='off';
    case 3
        BMW_Preferences.FitOptions.Display='detail';
end
switch get(handles.FitOptions_Algorithm,'Value')
    case 1
        BMW_Preferences.FitOptions.Algorithm='DE-MCMC';
    case 2
        BMW_Preferences.FitOptions.Algorithm='MH-MCMC';
    case 3
        BMW_Preferences.FitOptions.Algorithm='GA';
    case 4
        BMW_Preferences.FitOptions.Algorithm='SA';
    case 5
        BMW_Preferences.FitOptions.Algorithm='fmincon: sqp';
    case 6
        BMW_Preferences.FitOptions.Algorithm='fmincon: active-set';
    case 7
        BMW_Preferences.FitOptions.Algorithm='fmincon: interior-point';
    case 8
        BMW_Preferences.FitOptions.Algorithm='MADS';
    case 9
        BMW_Preferences.FitOptions.Algorithm='BADS';
end
% criteria
Criteria={};
if get(handles.LLH,'Value'), Criteria={Criteria,'LLH'}; end
if get(handles.AIC,'Value'), Criteria={Criteria,'AIC'}; end
if get(handles.AICc,'Value'), Criteria={Criteria,'AICc'}; end
if get(handles.BIC,'Value'), Criteria={Criteria,'BIC'}; end
if get(handles.DIC1,'Value'), Criteria={Criteria,'DIC1'}; end
if get(handles.DIC2,'Value'), Criteria={Criteria,'DIC2'}; end
if get(handles.DICs,'Value'), Criteria={Criteria,'DIC*'}; end
if get(handles.WAIC1,'Value'), Criteria={Criteria,'WAIC1'}; end
if get(handles.WAIC2,'Value'), Criteria={Criteria,'WAIC2'}; end
if get(handles.LME_BS,'Value'), Criteria={Criteria,'LME_BS'}; end
if get(handles.LME_GHM,'Value'), Criteria={Criteria,'LME_GHM'}; end
BMW_Preferences.Criteria=Criteria;
% model
ModelSpace={}; ModelVariants={};
if get(handles.Item_Limit,'Value'), ModelSpace={ModelSpace,'Item Limit'}; end
if get(handles.Standard_Mixture,'Value'), ModelSpace={ModelSpace,'Standard Mixture'}; end
if get(handles.Slots_plus_Averaging,'Value'), ModelSpace={ModelSpace,'Slots-plus-Averaging'}; end
if get(handles.Equal_Precision,'Value'), ModelSpace={ModelSpace,'Equal Precision'}; end
if get(handles.Variable_Precision,'Value'), ModelSpace={ModelSpace,'Variable Precision'}; end
if get(handles.Variable_Precision_with_Capacity,'Value'), ModelSpace={ModelSpace,'Variable Precision with Capacity'}; end
if get(handles.Category_Only,'Value'), ModelSpace={ModelSpace,'Category-Only'}; end
if get(handles.Category_Only_with_Capacity,'Value'), ModelSpace={ModelSpace,'Category-Only (with Capacity)'}; end
if get(handles.Fixed_Capacity_Central_Probe,'Value'), ModelSpace={ModelSpace,'Fixed Capacity (Central-Probe)'}; end
if get(handles.Fixed_Capacity_Full_Display,'Value'), ModelSpace={ModelSpace,'Fixed Capacity (Full-Display)'}; end
if get(handles.Fixed_Capacity_Single_Probe,'Value'), ModelSpace={ModelSpace,'Fixed Capacity (Single-Probe)'}; end
if get(handles.Signal_Detection,'Value'), ModelSpace={ModelSpace,'Signal Detection'}; end
if get(handles.Bias,'Value'), ModelVariants={ModelVariants,'Bias'}; end
if get(handles.BiasF, 'Value'), ModelVariants={ModelVariants, 'BiasF'}; end
if get(handles.PrecF,'Value'), ModelVariants={ModelVariants,'PrecF'}; end
if get(handles.Swap,'Value'), ModelVariants={ModelVariants,'Swap'}; end
if get(handles.ContinuousK,'Value'), ModelVariants={ModelVariants,'ContinuousK'}; end
if get(handles.ResponseNoise,'Value'), ModelVariants={ModelVariants,'ResponseNoise'}; end
if get(handles.Category_WI,'Value'), ModelVariants={ModelVariants,'Category (Within-Item)'}; end
if get(handles.Category_BI,'Value'), ModelVariants={ModelVariants,'Category (Between-Item)'}; end
if get(handles.Lapse,'Value'), ModelVariants={ModelVariants,'Lapse'}; end
if get(handles.Ensemble,'Value'), ModelVariants={ModelVariants,'Ensemble'}; end
BMW_Preferences.Model.Variants=ModelVariants;
BMW_Preferences.ModelSpace=ModelSpace;
