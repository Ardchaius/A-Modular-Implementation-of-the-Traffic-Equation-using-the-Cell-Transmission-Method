function varargout = CapandFlow(varargin)
% CAPANDFLOW MATLAB code for CapandFlow.fig
%      CAPANDFLOW, by itself, creates a new CAPANDFLOW or raises the existing
%      singleton*.
%
%      H = CAPANDFLOW returns the handle to a new CAPANDFLOW or the handle to
%      the existing singleton*.
%
%      CAPANDFLOW('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CAPANDFLOW.M with the given input arguments.
%
%      CAPANDFLOW('Property','Value',...) creates a new CAPANDFLOW or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before CapandFlow_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to CapandFlow_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help CapandFlow

% Last Modified by GUIDE v2.5 12-Jul-2019 23:06:44

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @CapandFlow_OpeningFcn, ...
                   'gui_OutputFcn',  @CapandFlow_OutputFcn, ...
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


% --- Executes just before CapandFlow is made visible.
function CapandFlow_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to CapandFlow (see VARARGIN)

% Choose default command line output for CapandFlow
handles.output = hObject;

if nargin<4 | ~isstruct(varargin{1})
else
    
    simSize = varargin{1}.Size; % Collect the input data
    road_Network = varargin{1}.Network;
    connections = varargin{2};
    segments = varargin{3};
    segLength = varargin{1}.segmentLength;
    capacity = varargin{4};
    flow = varargin{5};
    simWindow = varargin{6};
    dist = 0.7/(simSize); % Variables for the arrangment of the boxes
    boxSize = 0.7/(simSize);

    handles.SimSize = simSize; % Store variables globally
    handles.Network = road_Network; 
    handles.Connections = connections;
    handles.Segments = segments;
    handles.segmentLength = segLength;
    handles.SimWindow = simWindow;
    handles.Capacity = capacity;
    handles.Flow = flow;
    
    
    
    for i = 1:simSize
        for j = 1:simSize
            bID = i+(simSize-j)*simSize; % Button ID number
            tag = sprintf('pushbutton%i',bID); % Generate tag for button
            [loci, locj] = quorem(sym(bID-1), sym(simSize)); % Calculate location in standard matrix notation based on button number
            loci = loci+1;
            locj = locj+1;
            uicontrol('Parent',hObject,'Style','pushbutton','String',road_Network(loci,locj),'Units','normalized',...
                'Position',[0.15+dist*(i-1) 0.25+dist*(j-1) boxSize boxSize],'Visible','on',...
                'Callback',{@CaporFlow,hObject},'Tag', tag); % Create pushbutton
        end
    end
    
    uicontrol('Parent',hObject,'Style','pushbutton','String','Resume Simulation!',...
        'Units','normalized','Position',[0.3 0.05 0.4 0.15],'Visible','on',...
        'Callback',{@start_Simulation,hObject},'Tag', 'Start'); % Create start simulation button

end

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes CapandFlow wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = CapandFlow_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function CaporFlow(src,event,hObject)
% Function to ask user to edit Capacity or Flow of the chosen cell, or exit

handles = guidata(hObject);

CorF = questdlg('Do you want to edit the capacity or the flow of the selected segment?','Capacity or Flow','Capacity','Flow','Cancel','Cancel');

handles.CorF = CorF;

num = src.Tag;
num = str2double(num(11:end));
[loci, locj] = quorem(sym(num-1), sym(handles.SimSize)); % Calculate location in standard matrix notation based on button number
loci = loci+1;
locj = locj+1;

handles.Segi = loci;
handles.Segj = locj;

guidata(hObject,handles);

if ~strcmp(CorF,'Cancel')

    input = struct('Network',handles.Network,'simSize',handles.SimSize,...
        'segmentLength',handles.segmentLength,'CorF',handles.CorF,'Segi',handles.Segi,'Segj',handles.Segj);
    CapandFlow_Edit(input,handles.Connections,handles.Segments,handles.Capacity,handles.Flow,hObject);
    
end


function start_Simulation(src,event,hObject)
handles = guidata(hObject);

simWindow = handles.SimWindow;
simHandles = guidata(simWindow);
simHandles.Network = handles.Network;
simHandles.Connections = handles.Connections;
simHandles.Segments = handles.Segments;
simHandles.Capacity = handles.Capacity;
simHandles.Flow = handles.Flow;
guidata(simWindow,simHandles);
close;


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
simWindow = handles.SimWindow; % All buttons are enabled in the simulation window whenever the edit window is closed.
simHandles = guidata(simWindow);
set(simHandles.pushbutton1,'Enable','on');
set(simHandles.pushbutton2,'Enable','on');
set(simHandles.pushbutton3,'Enable','on');
set(simHandles.pushbutton4,'Enable','on');
set(simHandles.savebutton,'Enable','on');
set(simHandles.loadbutton,'Enable','on');
set(simHandles.pushbutton9,'Enable','on');
guidata(simWindow,simHandles);
% Hint: delete(hObject) closes the figure
delete(hObject);
