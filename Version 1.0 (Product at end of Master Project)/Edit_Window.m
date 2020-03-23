function varargout = Edit_Window(varargin)
% EDIT_WINDOW MATLAB code for Edit_Window.fig
%      EDIT_WINDOW, by itself, creates a new EDIT_WINDOW or raises the existing
%      singleton*.
%
%      H = EDIT_WINDOW returns the handle to a new EDIT_WINDOW or the handle to
%      the existing singleton*.
%
%      EDIT_WINDOW('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EDIT_WINDOW.M with the given input arguments.
%
%      EDIT_WINDOW('Property','Value',...) creates a new EDIT_WINDOW or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Edit_Window_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Edit_Window_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Edit_Window

% Last Modified by GUIDE v2.5 29-Mar-2019 09:56:23

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Edit_Window_OpeningFcn, ...
                   'gui_OutputFcn',  @Edit_Window_OutputFcn, ...
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


% --- Executes just before Edit_Window is made visible.
function Edit_Window_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Edit_Window (see VARARGIN)

% Choose default command line output for Edit_Window
handles.output = hObject;

if nargin<4 | ~isstruct(varargin{1})
else
    
    simSize = varargin{1}.Size; % Collect the input data
    road_Network = varargin{1}.Network;
    connections = varargin{2};
    segments = varargin{3};
    segLength = varargin{1}.segmentLength;
    inOut_T = varargin{4};
    simWindow = varargin{5};
    dist = 0.7/(simSize); % Variables for the arrangment of the boxes
    boxSize = 0.7/(simSize);

    handles.SimSize = simSize; % Store variables globally
    handles.Network = road_Network; 
    handles.Connections = connections;
    handles.Segments = segments;
    handles.segmentLength = segLength;
    handles.inOut_T = inOut_T;
    handles.SimWindow = simWindow;
    
    for i = 1:simSize
        for j = 1:simSize
            bID = i+(simSize-j)*simSize; % Button ID number
            tag = sprintf('pushbutton%i',bID); % Generate tag for button
            [loci, locj] = quorem(sym(bID-1), sym(simSize)); % Calculate location in standard matrix notation based on button number
            loci = loci+1;
            locj = locj+1;
            uicontrol('Parent',hObject,'Style','pushbutton','String',road_Network(loci,locj),'Units','normalized',...
                'Position',[0.15+dist*(i-1) 0.25+dist*(j-1) boxSize boxSize],'Visible','on',...
                'Callback',{@selection_screen,hObject},'Tag', tag); % Create pushbutton
        end
    end
    
    uicontrol('Parent',hObject,'Style','pushbutton','String','Resume Simulation!',...
        'Units','normalized','Position',[0.3 0.05 0.4 0.15],'Visible','on',...
        'Callback',{@start_Simulation,hObject},'Tag', 'Start'); % Create start simulation button

end

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Edit_Window wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Edit_Window_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function selection_screen(src,event,hObject)

handles = guidata(hObject);
simSize = handles.SimSize; % Get global varibales
road_Network = handles.Network;
connections = handles.Connections;
segments = handles.Segments;
segLength = handles.segmentLength;
inOut_T = handles.inOut_T;

slr = sprintf('Straight Left-Right %c',char(9472)); % Straight road section Left-Right
stb = sprintf('Straight Top-Bottom %c',char(9474)); % Straight road section Top-Bottom
clb = sprintf('Corner Left-Bottom %c',char(9488)); % Corner road section, left to bottom
cbr = sprintf('Corner Bottom-Right %c', char(9484)); % Corner road section, bottom to right
crt = sprintf('Corner Right-Top %c', char(9492)); % Corner road section, right to top
ctl = sprintf('Corner Top-Left %c', char(9496)); % Corner road section, top to left
isl = sprintf('Intersection from Left %c',char(9508)); % Intersection from Left
isb = sprintf('Intersection from Bottom %c',char(9516)); % Intersection from Bottom
isr = sprintf('Intersection from Right %c',char(9500)); % Intersection from Right
ist = sprintf('Intersection from Top %c',char(9524)); % Intersection from Top
rda = sprintf('Roundabout %c',char(9711)); % Roundabout with entrances and exits in all directions
tlc = sprintf('Traffic Light Crossing %c', char(9532)); % Traffic light 4-way crossing
empty = 'Empty cell';
entex = 'A spawning and exiting point for cars EX';
list = {slr, stb, clb, cbr, crt, ctl, isl, isb, isr, ist, rda, tlc, empty, entex};
[index,~] = listdlg('PromptString','Please Sselect a road type for the road section:',...
    'SelectionMode','single','ListString',list); % Pop up dialog list box
if isempty(index) % Line of code to avoid annoying ding sound from matlab when pressing cancel, switch needs a scalar
    index = 0;
end

num = src.Tag;
num = str2double(num(11:end));
[loci, locj] = quorem(sym(num-1), sym(simSize)); % Calculate location in standard matrix notation based on button number
loci = loci+1;
locj = locj+1;

switch index
    
    case 1 % Left - Right
        src.String = char(9472);
        connections(loci,locj) = {[0,1,0,1]};
        segments(loci,locj) = {zeros(segLength,2)}; % Straight road sections stored as: [ml,mr] with both being vertical vectors
        inOut_T(loci,locj) = {zeros(2,4)};
    case 2 % Top - Bottom
        src.String = char(9474);
        connections(loci,locj) = {[1,0,1,0]};
        segments(loci,locj) = {zeros(segLength,2)}; % Straight road sections stored as: [ml,mr] with both being vertical vectors
        inOut_T(loci,locj) = {zeros(2,4)};
    case 3 % Left - Bottom
        src.String = char(9488);
        connections(loci,locj) = {[0,0,1,1]};
        segments(loci,locj) = {zeros(segLength-2,2)}; % Corner road sections stored as: [ml,mr] with both being vertical vectors and 2 shorter than a straight road due to curve
        inOut_T(loci,locj) = {zeros(2,4)};
    case 4 % Bottom - Right
        src.String = char(9484);
        connections(loci,locj) = {[0,1,1,0]};
        segments(loci,locj) = {zeros(segLength-2,2)}; % Corner road sections stored as: [ml,mr] with both being vertical vectors and 2 shorter than a straight road due to curve
        inOut_T(loci,locj) = {zeros(2,4)};
    case 5 % Right - Top
        src.String = char(9492);
        connections(loci,locj) = {[1,1,0,0]};
        segments(loci,locj) = {zeros(segLength-2,2)}; % Corner road sections stored as: [ml,mr] with both being vertical vectors and 2 shorter than a straight road due to curve
        inOut_T(loci,locj) = {zeros(2,4)};
    case 6 % Top - Left
        src.String = char(9496);
        connections(loci,locj) = {[1,0,0,1]};
        segments(loci,locj) = {zeros(segLength-2,2)}; % Corner road sections stored as: [ml,mr] with both being vertical vectors and 2 shorter than a straight road due to curve
        inOut_T(loci,locj) = {zeros(2,4)};
    case 7 % Intersection Left
        src.String = char(9508);
        connections(loci,locj) = {[1,0,1,1]};
        segments(loci,locj) = {zeros(segLength+segLength/2-1/2,2)}; % Corner road sections stored as: [ml,mr] with both being vertical vectors and 2 shorter than a straight road due to curve
        inOut_T(loci,locj) = {zeros(2,4)};
    case 8 % Intersection Bottom
        src.String = char(9516);
        connections(loci,locj) = {[0,1,1,1]};
        segments(loci,locj) = {zeros(segLength+segLength/2-1/2,2)}; % Corner road sections stored as: [ml,mr] with both being vertical vectors and 2 shorter than a straight road due to curve
        inOut_T(loci,locj) = {zeros(2,4)};
    case 9 % Intersection Right
        src.String = char(9500);
        connections(loci,locj) = {[1,1,1,0]};
        segments(loci,locj) = {zeros(segLength+segLength/2-1/2,2)}; % Corner road sections stored as: [ml,mr] with both being vertical vectors and 2 shorter than a straight road due to curve
        inOut_T(loci,locj) = {zeros(2,4)};
    case 10 % Intersection Top
        src.String = char(9524);
        connections(loci,locj) = {[1,1,0,1]};
        segments(loci,locj) = {zeros(segLength+segLength/2-1/2,2)}; % Corner road sections stored as: [ml,mr] with both being vertical vectors and 2 shorter than a straight road due to curve
        inOut_T(loci,locj) = {zeros(2,4)};
    case 11 % Roundabout
        src.String = char(9711);
        connections(loci,locj) = {[1,1,1,1]};
        segments(loci,locj) = {zeros(segLength-1,4)}; % Corner road sections stored as: [ml,mr] with both being vertical vectors and 2 shorter than a straight road due to curve
        inOut_T(loci,locj) = {zeros(2,4)};
    case 12 % 4-Way Traffic Light
        src.String = char(9532);
        connections(loci,locj) = {[1,1,1,1]};
        segments(loci,locj) = {[zeros(segLength-1,4);1,0,0,0]}; % Corner road sections stored as: [ml,mr] with both being vertical vectors and 2 shorter than a straight road due to curve
        inOut_T(loci,locj) = {zeros(2,4)};
    case 13
        src.String = '';
        connections(loci,locj) = {[0,0,0,0]};
        segments(loci,locj) = {0};
        inOut_T(loci,locj) = {[]};
    case 14 % Input
        src.String = 'EX';
        connections(loci,locj) = {[1,1,1,1]};
        segments(loci,locj) = {0};
        inOut_T(loci,locj) = {zeros(2,4)};
    otherwise
    
end

road_Network(loci,locj) = src.String;
handles.Network = road_Network; % Update global network
handles.Connections = connections;
handles.Segments = segments;
handles.inOut_T = inOut_T;

guidata(hObject,handles); % Update global handles

function start_Simulation(src,event,hObject)
handles = guidata(hObject);
valid = check_Validity(handles.Network,handles.Connections,handles.SimSize); % Check validity of edited system
if valid % Update the global variables in the simulation window and clsoe this window
    simWindow = handles.SimWindow;
    simHandles = guidata(simWindow);
    simHandles.Network = handles.Network;
    simHandles.Connections = handles.Connections;
    simHandles.Segments = handles.Segments;
    simHandles.InOut_T = handles.inOut_T;
    guidata(simWindow,simHandles);
    close;
end


function valid = check_Validity(road_Network,connections,simSize)

invalid_Pieces = zeros(0);
% neighbour_connections, Matrix to indicate what connections the neigbouring segment has. row 1 is above, 2- to the right, 3- below, 4- to the left

for i = 1:simSize
    for j = 1:simSize
        
        segment = road_Network(i,j);
        
        if strcmp(segment,'') || strcmp(segment,'EX')
            continue;
        end
        
        local_connections = connections{i,j};
        
        if (i == 1) && (local_connections(1) == 1) % Check whether the segment is in the top section and has connections further upwards
            invalid_Pieces(end+1,:) = [i,j,1];
        elseif (local_connections(1) == 1)
            neighbour_connections = connections{i-1,j};
            if (neighbour_connections(3) == 0)
                invalid_Pieces(end+1,:) = [i,j,1];
            end
        end
        
        if (i == simSize) && (local_connections(3) == 1) % Check whether the segment is in the bottom section and has connections further downwards
            invalid_Pieces(end+1,:) = [i,j,3];
        elseif (local_connections(3) == 1)
            neighbour_connections = connections{i+1,j};
            if (neighbour_connections(1) == 0)
                invalid_Pieces(end+1,:) = [i,j,3];
            end
        end
        
        if (j == 1) && (local_connections(4) == 1)
            invalid_Pieces(end+1,:) = [i,j,4];
        elseif (local_connections(4) == 1)
            neighbour_connections = connections{i,j-1};
            if (neighbour_connections(2) == 0)
                invalid_Pieces(end+1,:) = [i,j,4];
            end
        end
        
        if (j == simSize) && (local_connections(2) == 1)
            invalid_Pieces(end+1,:) = [i,j,2];
        elseif (local_connections(2) == 1)
            neighbour_connections = connections{i,j+1};
            if (neighbour_connections(4) == 0)
                invalid_Pieces(end+1,:) = [i,j,2];
            end
        end
        
    end
end

if isempty(invalid_Pieces)
        valid = 1;
    return;
end
        
valid = 0;
error_msg = 'The chosen system is not valid. Invalid pieces are listed below:';

[invalid_num,~] = size(invalid_Pieces);

for i = 1:invalid_num
    switch invalid_Pieces(i,3)
        case 1
            temp1 = 'above.';
        case 2
            temp1 = 'to the right.';
        case 3
            temp1 = 'below.';
        case 4
            temp1 = 'to the left.';
    end
    temp2 = sprintf('The segment in position (%i,%i) has no valid connection ',invalid_Pieces(i,1),invalid_Pieces(i,2));
    error_msg = [error_msg newline temp2 temp1];
end

errordlg(error_msg);


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
guidata(simWindow,simHandles);
% Hint: delete(hObject) closes the figure
delete(hObject);
