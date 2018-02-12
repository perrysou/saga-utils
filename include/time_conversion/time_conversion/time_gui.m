function varargout = time_gui(varargin)
% TIME_GUI M-file for time_gui.fig
%      TIME_GUI, by itself, creates a new TIME_GUI or raises the existing
%      singleton*.
%
%      H = TIME_GUI returns the handle to a new TIME_GUI or the handle to
%      the existing singleton*.
%
%      TIME_GUI('Property','Value',...) creates a new TIME_GUI using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to time_gui_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      TIME_GUI('CALLBACK') and TIME_GUI('CALLBACK',hObject,...) call the
%      local function named CALLBACK in TIME_GUI.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help time_gui

% Last Modified by GUIDE v2.5 19-Jun-2003 14:18:43

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @time_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @time_gui_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
                   'gui_Callback',   []);
if nargin & isstr(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before time_gui is made visible.
function time_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for time_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes time_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = time_gui_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in leap_seconds.
function leap_seconds_Callback(hObject, eventdata, handles)
% hObject    handle to leap_seconds (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of leap_seconds


% --- Executes when figure1 window is resized.
function figure1_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


