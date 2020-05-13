function varargout = BMW_About(varargin)

% Last Modified by GUIDE v2.5 17-Apr-2020 13:31:14

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @BMW_About_OpeningFcn, ...
                   'gui_OutputFcn',  @BMW_About_OutputFcn, ...
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


% --- Executes just before BMW_About is made visible.
function BMW_About_OpeningFcn(hObject, eventdata, handles, varargin)

handles.output = hObject;
set(handles.axes1,'visible','off')
axes(handles.axes1);
image1=imread('BMW_icon.png');
imshow(image1);
set(handles.axes2,'visible','off')
axes(handles.axes2);
image2=imread('MAC_HotPot.jpg');
imshow(image2);
% Update handles structure
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = BMW_About_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
