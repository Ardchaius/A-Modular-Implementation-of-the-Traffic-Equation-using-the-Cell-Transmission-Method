function varargout = Simulation_Window(varargin)
% SIMULATION_WINDOW MATLAB code for Simulation_Window.fig
%      SIMULATION_WINDOW, by itself, creates a new SIMULATION_WINDOW or raises the existing
%      singleton*.
%
%      H = SIMULATION_WINDOW returns the handle to a new SIMULATION_WINDOW or the handle to
%      the existing singleton*.
%
%      SIMULATION_WINDOW('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SIMULATION_WINDOW.M with the given input arguments.
%
%      SIMULATION_WINDOW('Property','Value',...) creates a new SIMULATION_WINDOW or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Simulation_Window_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Simulation_Window_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Simulation_Window

% Last Modified by GUIDE v2.5 12-Jul-2019 23:04:21

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Simulation_Window_OpeningFcn, ...
                   'gui_OutputFcn',  @Simulation_Window_OutputFcn, ...
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


% --- Executes just before Simulation_Window is made visible.
function Simulation_Window_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Simulation_Window (see VARARGIN)

% Choose default command line output for Simulation_Window
handles.output = hObject;

if nargin<4 | ~isstruct(varargin{1})
else
    simSize = varargin{1}.Size; % Collect the size
    dist = 0.7/(simSize); % Variables for the arrangment of the boxes
    boxSize = 0.7/(simSize);
    
    road_Network = strings(simSize); % Create an empty array to store network for calculation
    road_Network(:) = ''; % Initialize as empty segments

    handles.SimSize = simSize; % Store size globally
    handles.Network = road_Network; % Store network globally
    handles.Connections = cell(simSize); % Store connections globally
    connections = handles.Connections;
    
    for i = 1:simSize
        for j = 1:simSize
            tag = sprintf('pushbutton%i',i+(simSize-j)*simSize); % Generate tag for button
            uicontrol('Parent',hObject,'Style','pushbutton','String','','Units','normalized',...
                'Position',[0.15+dist*(i-1) 0.25+dist*(j-1) boxSize boxSize],'Visible','on',...
                'Callback',{@selection_screen,hObject},'Tag', tag); % Create pushbutton in appropriate location
            connections(i,j) = {[0,0,0,0]};
        end
    end
    
    handles.Connections = connections;
    
    uicontrol('Parent',hObject,'Style','pushbutton','String','Start Simulation!',...
        'Units','normalized','Position',[0.3 0.05 0.4 0.15],'Visible','on',...
        'Callback',{@start_Simulation,hObject},'Tag', 'Start'); % Create start simulation button

end

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Simulation_Window wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Simulation_Window_OutputFcn(hObject, eventdata, handles) 
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
list = {slr, stb, clb, cbr, crt, ctl, isl, isb, isr, ist, rda, tlc, empty, entex}; % Create list for dialog box
[index,~] = listdlg('PromptString','Please Select a road type for the road section:',...
    'SelectionMode','single','ListString',list); % Pop up dialog list box
if isempty(index) % Line of code to avoid annoying ding sound from matlab when pressing cancel, switch needs a scalar
    index = 0;
end

num = src.Tag; % Get the number of the button that is being pushed
num = str2double(num(11:end));
[loci, locj] = quorem(sym(num-1), sym(simSize)); % Calculate location in standard matrix notation based on button number
loci = loci+1;
locj = locj+1;

switch index
    
    case 1 % Left - Right
        src.String = char(9472); % Change the text on the butoton that is being edited
        connections(loci,locj) = {[0,1,0,1]}; % Add appropriate connections, 1 is north, 2 is east, 3 is south and 4 is west
        % 1 indicates the segment has a road coming from or going to that
        % particular direction.
    case 2 % Top - Bottom
        src.String = char(9474);
        connections(loci,locj) = {[1,0,1,0]};
    case 3 % Left - Bottom
        src.String = char(9488);
        connections(loci,locj) = {[0,0,1,1]};
    case 4 % Bottom - Right
        src.String = char(9484);
        connections(loci,locj) = {[0,1,1,0]};
    case 5 % Right - Top
        src.String = char(9492);
        connections(loci,locj) = {[1,1,0,0]};
    case 6 % Top - Left
        src.String = char(9496);
        connections(loci,locj) = {[1,0,0,1]};
    case 7 % Intersection Left
        src.String = char(9508);
        connections(loci,locj) = {[1,0,1,1]};
    case 8 % Intersection Bottom
        src.String = char(9516);
        connections(loci,locj) = {[0,1,1,1]};
    case 9 % Intersection Right
        src.String = char(9500);
        connections(loci,locj) = {[1,1,1,0]};
    case 10 % Intersection Top
        src.String = char(9524);
        connections(loci,locj) = {[1,1,0,1]};
    case 11 % Roundabout
        src.String = char(9711);
        connections(loci,locj) = {[1,1,1,1]};
    case 12 % 4-Way Traffic Light
        src.String = char(9532);
        connections(loci,locj) = {[1,1,1,1]};
    case 13
        src.String = '';
        connections(loci,locj) = {[0,0,0,0]};
    case 14 % Input
        src.String = 'EX';
        connections(loci,locj) = {[1,1,1,1]};
    otherwise 
    
end

road_Network(loci,locj) = src.String;
handles.Network = road_Network; % Update global network
handles.Connections = connections;

guidata(hObject,handles); % Update global handles

function start_Simulation(src,event,hObject)
handles = guidata(hObject);
valid = check_Validity(handles.Network,handles.Connections,handles.SimSize); % Check the validity of the system
if valid
    input = struct('Network',handles.Network,'Size',handles.SimSize);
    close;
    Simulation_Plot(input,handles.Connections); % Start simulation window
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
        elseif (local_connections(1) == 1) % Check if the segment to the north is a valid connection
            neighbour_connections = connections{i-1,j};
            if (neighbour_connections(3) == 0)
                invalid_Pieces(end+1,:) = [i,j,1];
            end
        end
        
        if (i == simSize) && (local_connections(3) == 1) % Check whether the segment is in the bottom section and has connections further downwards
            invalid_Pieces(end+1,:) = [i,j,3];
        elseif (local_connections(3) == 1) % Check if the segment to the south is a valid connection
            neighbour_connections = connections{i+1,j};
            if (neighbour_connections(1) == 0)
                invalid_Pieces(end+1,:) = [i,j,3];
            end
        end
        
        if (j == 1) && (local_connections(4) == 1) % Check whether the segment is in the left section and has connections further to the left
            invalid_Pieces(end+1,:) = [i,j,4];
        elseif (local_connections(4) == 1) % Check if the segment to the left is a valid connection
            neighbour_connections = connections{i,j-1};
            if (neighbour_connections(2) == 0)
                invalid_Pieces(end+1,:) = [i,j,4];
            end
        end
        
        if (j == simSize) && (local_connections(2) == 1) % Check whether the segment is in the right section and has connections further to the right
            invalid_Pieces(end+1,:) = [i,j,2];
        elseif (local_connections(2) == 1) % Check if the segment to the right is a valid connection
            neighbour_connections = connections{i,j+1};
            if (neighbour_connections(4) == 0)
                invalid_Pieces(end+1,:) = [i,j,2];
            end
        end
        
    end
end

if isempty(invalid_Pieces) % Check if the system has invalid pieces, return otherwise
        valid = 1;
    return;
end
        
valid = 0;
error_msg = 'The chosen system is not valid. Invalid pieces are listed below:';

[invalid_num,~] = size(invalid_Pieces);

for i = 1:invalid_num % Create appropriate error message for each invalid piece
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

errordlg(error_msg); % Display error message
