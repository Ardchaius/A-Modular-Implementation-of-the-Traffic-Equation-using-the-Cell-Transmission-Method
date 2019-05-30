function varargout = Start_Up_Window(varargin)
% START_UP_WINDOW MATLAB code for Start_Up_Window.fig
%      START_UP_WINDOW, by itself, creates a new START_UP_WINDOW or raises the existing
%      singleton*.
%
%      H = START_UP_WINDOW returns the handle to a new START_UP_WINDOW or the handle to
%      the existing singleton*.
%
%      START_UP_WINDOW('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in START_UP_WINDOW.M with the given input arguments.
%
%      START_UP_WINDOW('Property','Value',...) creates a new START_UP_WINDOW or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Start_Up_Window_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Start_Up_Window_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Start_Up_Window

% Last Modified by GUIDE v2.5 14-Mar-2019 23:09:55

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Start_Up_Window_OpeningFcn, ...
                   'gui_OutputFcn',  @Start_Up_Window_OutputFcn, ...
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


% --- Executes just before Start_Up_Window is made visible.
function Start_Up_Window_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Start_Up_Window (see VARARGIN)

% Choose default command line output for Start_Up_Window
handles.output = hObject;

set(handles.text2, 'String', 'Please enter the size of the system you wish to simulate:');

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Start_Up_Window wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Start_Up_Window_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



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


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

simSize = str2double(get(handles.edit1, 'String'));

if (1 <= simSize) && (simSize <= 10)
    input = struct('Size', simSize);
    close;
    Simulation_Window(input);
else
    errordlg('The size of the system you have chosen is not suported, please ensure to choose a square size between 1 and 10.','Invalid size of System.','modal');
end
