function varargout = Simulation_Plot(varargin)
% SIMULATION_PLOT MATLAB code for Simulation_Plot.fig
%      SIMULATION_PLOT, by itself, creates a new SIMULATION_PLOT or raises the existing
%      singleton*.
%
%      H = SIMULATION_PLOT returns the handle to a new SIMULATION_PLOT or the handle to
%      the existing singleton*.
%
%      SIMULATION_PLOT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SIMULATION_PLOT.M with the given input arguments.
%
%      SIMULATION_PLOT('Property','Value',...) creates a new SIMULATION_PLOT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Simulation_Plot_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Simulation_Plot_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Simulation_Plot

% Last Modified by GUIDE v2.5 12-Jul-2019 21:28:33

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Simulation_Plot_OpeningFcn, ...
                   'gui_OutputFcn',  @Simulation_Plot_OutputFcn, ...
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


% --- Executes just before Simulation_Plot is made visible.
function Simulation_Plot_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Simulation_Plot (see VARARGIN)

% Choose default command line output for Simulation_Plot
handles.output = hObject;

if nargin<4 | ~isstruct(varargin{1})
else
    
    tempstring = sprintf('Auto Play %c',char(9654)); % Intialize autoplay button
    set(handles.pushbutton2, 'String', tempstring);
    
    road_Network = varargin{1}.Network; % Retrieve the Network and Size from input
    simSize = varargin{1}.Size;
    connections = varargin{2};
    
    handles.InitCap = 3; % Initial capacity
    handles.InitFlow = 3; % Initial flow
    handles.Network = road_Network; % Save input in gui handles and initialize Segment cell
    handles.Size = simSize;
    handles.Connections = connections;
    set(0,'UserData',0);
    handles.Segments = cell(simSize); % For time t
    handles.InOut_T = cell(simSize); % Input and output for each segment in each direction
    % Each entry in the cell consists of a 2x4 matrix. Row one is the
    % output in each direction and row 2 is the input in each direction. 
    % The column indicates the direction as
    % follows: column 1 = top, column 2 = right, column 3 = bottom, column 4 = left
    % For time t
    handles.Capacity = cell(simSize); % Capacity for each cell in a given segment
    handles.Flow = cell(simSize); % Maximum flow rate for each cell in a given segment
    handles.InOut_T1 = cell(simSize);
    handles.InOut_T1(:) = {zeros(2,4)};
    handles.InOut_C = cell(simSize);
    handles.InOut_C(:) = {handles.InitCap*ones(2,4)};
    handles.InOut_F = cell(simSize);
    handles.InOut_F(:) = {handles.InitFlow*ones(2,4)};
    handles.SegmentLength = 9; % Static for now, later controlled by user, has to be odd and equal t or greater than 7.
    handles.InputVolume = 'random'; % As with Segment length, will later be controlled by user 'random' randomizes input
    axes(handles.axes1); % Store the axes for plotting of the simulation
    
    handles.cMap = [1 1 1;
        (0:0.001:1)',(1:-0.001:0)',zeros(1001,1)]; % Custom colour map for simulation
    
    guidata(hObject,handles); % Save global variables
    
    initializeSegments(hObject); % Initialize the Segments cell
    handles = guidata(hObject);
    
end

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Simulation_Plot wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Simulation_Plot_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
% Advance simulation one time step and plot the data in the Segments cell array
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Simulation_Single_Time_Step(hObject);



% --- Executes on button press in pushbutton2.
% Autoplay simulation
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(0,'UserData',~get(0,'UserData'));
if get(0,'UserData') % Is true when auto simulation is started, disables all other buttons
    tempstring = sprintf('Auto Play %c',char(9208)); % Change text on autoplay button to have pause symbol
    set(handles.pushbutton4,'Enable','off');
    set(handles.pushbutton3,'Enable','off');
    set(handles.pushbutton1,'Enable','off');
    set(handles.savebutton,'Enable','off');
    set(handles.loadbutton,'Enable','off');
    set(handles.pushbutton9,'Enable','off');
else % Is true when auto zimulation is stopped, enables all other buttons
    tempstring = sprintf('Auto Play %c',char(9654)); % Change text on autoplay button to have play symbol
    set(handles.pushbutton4,'Enable','on');
    set(handles.pushbutton3,'Enable','on');
    set(handles.pushbutton1,'Enable','on');
    set(handles.savebutton,'Enable','on');
    set(handles.loadbutton,'Enable','on');
    set(handles.pushbutton9,'Enable','on');
end

set(handles.pushbutton2, 'String', tempstring);

while get(0,'UserData')
    Simulation_Single_Time_Step(hObject);
end


% --- Executes on button press in pushbutton3.
% Edit the road network and resume the simulation
% All exchanged pieces start empty
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Disable all buttons in the Simulation Window
set(handles.pushbutton1,'Enable','off');
set(handles.pushbutton2,'Enable','off');
set(handles.pushbutton3,'Enable','off');
set(handles.pushbutton4,'Enable','off');
set(handles.savebutton,'Enable','off');
set(handles.loadbutton,'Enable','off');
set(handles.pushbutton9,'Enable','off');
guidata(hObject,handles); % Update global variables
% Format input for edit window
input = struct('Network',handles.Network,'Size',handles.Size,'segmentLength',handles.SegmentLength);
Edit_Window(input,handles.Connections,handles.Segments,handles.InOut_T,hObject); % Open edit window



% --- Executes on button press in pushbutton4. Reset the simulation
% Reset the simulation back to a zero state
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

initializeSegments(hObject); % Initialize the Segments cell, in this case resets to null
handles = guidata(hObject);
guidata(hObject,handles); % Update global variables



% --- Executes on button press in savebutton.
% Save the simulation to text file
function savebutton_Callback(hObject, eventdata, handles)
% hObject    handle to savebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[handles.FileName, ~] = uiputfile('*.txt'); % Open dialog to choose save file
guidata(hObject,handles); % Update global variables

save_Network(handles); % Save function



% --- Executes on button press in loadbutton.
% Load a previously saved simulation
function loadbutton_Callback(hObject, eventdata, handles)
% hObject    handle to loadbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.FileName = uigetfile('*.txt'); % Open dialog to choose file to load
guidata(hObject,handles); % Update global variables

load_Network(hObject,handles); % Load function


% --- Executes on button press in pushbutton9.
% Edit cell capacity and flow rate
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Disable all buttons in the Simulation Window
set(handles.pushbutton1,'Enable','off');
set(handles.pushbutton2,'Enable','off');
set(handles.pushbutton3,'Enable','off');
set(handles.pushbutton4,'Enable','off');
set(handles.savebutton,'Enable','off');
set(handles.loadbutton,'Enable','off');
set(handles.pushbutton9,'Enable','off');
guidata(hObject,handles); % Update global variables
% Format input for edit window
input = struct('Network',handles.Network,'Size',handles.Size,'segmentLength',handles.SegmentLength);
CapandFlow(input,handles.Connections,handles.Segments,handles.Capacity,handles.Flow,hObject); % Open edit window



function initializeSegments(hObject)

% Retrieve global variables
handles = guidata(hObject);
simSize = handles.Size;
segments = handles.Segments;
inOut_T = handles.InOut_T;
road_Network = handles.Network;
segLength = handles.SegmentLength;
IV = handles.InputVolume;
capacity = handles.Capacity;
flow = handles.Flow;
initCap = handles.InitCap;
initFlow = handles.InitFlow;

for i = 1:simSize % Loop for each road segment in the simulation and assign appropriate variables
    for j = 1:simSize
        
        switch road_Network(i,j)
           
            case ''
                segments(i,j) = {0};
            case 'EX'
                if strcmp(IV,'random')
                    tempIV = zeros(1,4);
                    for k = 1:4
                        tempIV(k) = randi([0,3]);
                    end
                    segments(i,j) = {0};
                    inOut_T(i,j) = {[tempIV;zeros(1,4)]};
                else
                    segments(i,j) = {0};
                    inOut_T(i,j) = {[IV,IV,IV,IV;0,0,0,0]};
                end
                capacity(i,j) = {inf*ones(2,4)};
                flow(i,j) = {initFlow*ones(2,4)};
            case char(9472) % L-R
                segments(i,j) = {zeros(segLength,2)}; % Straight road sections stored as: [ml,mr] with both being vertical vectors
                inOut_T(i,j) = {zeros(2,4)};
                capacity(i,j) = {initCap*ones(segLength,2)};
                flow(i,j) = {initFlow*ones(segLength,2)};
            case char(9474) % T-B
                segments(i,j) = {zeros(segLength,2)}; % Straight road sections stored as: [ml,mr] with both being vertical vectors
                inOut_T(i,j) = {zeros(2,4)};
                capacity(i,j) = {initCap*ones(segLength,2)};
                flow(i,j) = {initFlow*ones(segLength,2)};
            case char(9488) % L-B
                segments(i,j) = {zeros(segLength-2,2)}; % Corner road sections stored as: [ml,mr] with both being vertical vectors and 2 shorter than a straight road due to curve
                inOut_T(i,j) = {zeros(2,4)};
                capacity(i,j) = {initCap*ones(segLength-2,2)};
                flow(i,j) = {initFlow*ones(segLength-2,2)};
            case char(9484) % B-R
                segments(i,j) = {zeros(segLength-2,2)}; % Corner road sections stored as: [ml,mr] with both being vertical vectors and 2 shorter than a straight road due to curve
                inOut_T(i,j) = {zeros(2,4)};
                capacity(i,j) = {initCap*ones(segLength-2,2)};
                flow(i,j) = {initFlow*ones(segLength-2,2)};
            case char(9492) % R-T
                segments(i,j) = {zeros(segLength-2,2)}; % Corner road sections stored as: [ml,mr] with both being vertical vectors and 2 shorter than a straight road due to curve
                inOut_T(i,j) = {zeros(2,4)};
                capacity(i,j) = {initCap*ones(segLength-2,2)};
                flow(i,j) = {initFlow*ones(segLength-2,2)};
            case char(9496) % T-L
                segments(i,j) = {zeros(segLength-2,2)}; % Corner road sections stored as: [ml,mr] with both being vertical vectors and 2 shorter than a straight road due to curve
                inOut_T(i,j) = {zeros(2,4)};
                capacity(i,j) = {initCap*ones(segLength-2,2)};
                flow(i,j) = {initFlow*ones(segLength-2,2)};
            case char(9508) % Intersection from Left
                segments(i,j) = {zeros(segLength+segLength/2-1/2,2)}; % Intersections are stored in the following fashion: [ml,mr;sl,sr;t2,t1]
                inOut_T(i,j) = {zeros(2,4)};
                capacity(i,j) = {initCap*ones(segLength+segLength/2-1/2,2)};
                flow(i,j) = {initFlow*ones(segLength+segLength/2-1/2,2)};
            case char(9516) % Intersection from Bottom
                segments(i,j) = {zeros(segLength+segLength/2-1/2,2)}; % Intersections are stored in the following fashion: [ml,mr;sl,sr;t2,t1]
                inOut_T(i,j) = {zeros(2,4)};
                capacity(i,j) = {initCap*ones(segLength+segLength/2-1/2,2)};
                flow(i,j) = {initFlow*ones(segLength+segLength/2-1/2,2)};
            case char(9500) % Intersection from Right
                segments(i,j) = {zeros(segLength+segLength/2-1/2,2)}; % Intersections are stored in the following fashion: [ml,mr;sl,sr;t2,t1]
                inOut_T(i,j) = {zeros(2,4)};
                capacity(i,j) = {initCap*ones(segLength+segLength/2-1/2,2)};
                flow(i,j) = {initFlow*ones(segLength+segLength/2-1/2,2)};
            case char(9524) % Intersection from Top
                segments(i,j) = {zeros(segLength+segLength/2-1/2,2)}; % Intersections are stored in the following fashion: [ml,mr;sl,sr;t2,t1]
                inOut_T(i,j) = {zeros(2,4)};
                capacity(i,j) = {initCap*ones(segLength+segLength/2-1/2,2)};
                flow(i,j) = {initFlow*ones(segLength+segLength/2-1/2,2)};
            case char(9711) % Roundabout
                segments(i,j) = {zeros(segLength-1,4)}; % Roundabouts stored as [mt,mr,mb,ml;stl,srl,sbl,sll;str,srr,sbr,slr];
                inOut_T(i,j) = {zeros(2,4)};
                capacity(i,j) = {initCap*ones(segLength-1,4)};
                flow(i,j) = {initFlow*ones(segLength-1,4)};
            case char(9532) % Traffic light
                segments(i,j) = {[zeros(segLength-1,4);1,0,0,0]}; % 4-Way Traffic light intersection saved as [tl,rl,bl,ll;tr,rr,br,lr;tc,rc,bc,lc;timer,0,0,0] Here the right lane includes the cell in the intersection itself
                inOut_T(i,j) = {zeros(2,4)};
                capacity(i,j) = {initCap*ones(segLength-1,4)};
                flow(i,j) = {initFlow*ones(segLength-1,4)};
                
        end
        
    end
end
handles.Segments = segments;
handles.InOut_T = inOut_T;
handles.InputVolume = IV;
handles.Capacity = capacity;
handles.Flow = flow;

visualize_Network(handles); % Visualize the initial state

guidata(hObject, handles); % Update global variables


function Simulation_Single_Time_Step(hObject)

handles = guidata(hObject); % Retrieve global variables

simSize = handles.Size;
segments = handles.Segments;
road_Network = handles.Network;
length = handles.SegmentLength;
inOut_T = handles.InOut_T;
inOut_T1 = handles.InOut_T1;
IV = handles.InputVolume; % Input volume for input cells
capacity = handles.Capacity;
flow = handles.Flow;
inOut_C = handles.InOut_C;
inOut_F = handles.InOut_F;

for i = 1:simSize
    for j = 1:simSize
        
        temp = segments{i,j}; % Saving data matrix to temporary variable for alteration before passing to functions
        temp_C = capacity{i,j};
        temp_F = flow{i,j};
        
        switch road_Network(i,j)
            
            case 'EX'
                if strcmp(IV,'random')
                    tempIV = zeros(1,4);
                    for k = 1:4
                        tempIV(k) = randi([0,3]);
                    end
                    segments(i,j) = {0};
                    inOut_T1(i,j) = {[tempIV;zeros(1,4)]};
                else
                    segments(i,j) = {0};
                    inOut_T1(i,j) = {[IV,IV,IV,IV;0,0,0,0]};
                end
            case char(9472) % Straight L-R
                temp_left = inOut_T{i,j-1}; % Format input data
                temp_right = inOut_T{i,j+1};
                temp_left_C = inOut_C{i,j-1};
                temp_right_C = inOut_C{i,j+1};
                temp_left_F = inOut_F{i,j-1};
                temp_right_F = inOut_F{i,j+1};
                mlt = [temp_right(1,4);temp(:,1); temp_left(2,2)];
                mrt = [temp_left(1,2);temp(:,2);temp_right(2,4)];
                mlt_C = [temp_right_C(1,4);temp_C(:,1); temp_left_C(2,2)];
                mrt_C = [temp_left_C(1,2);temp_C(:,2);temp_right_C(2,4)];
                mlt_F = [temp_right_F(1,4);temp_F(:,1); temp_left_F(2,2)];
                mrt_F = [temp_left_F(1,2);temp_F(:,2);temp_right_F(2,4)];
                [mlt1,mrt1] = GUI_Straight_CF(mlt,mrt,mlt_C,mrt_C,mlt_F,mrt_F,length); % Calculation function
                inOut_T1(i,j) = {[0,mrt1(end),0,mlt1(end);... % Format output data
                    0,mlt1(1),0,mrt1(1)]};
                segments(i,j) = {[mlt1,mrt1]};
            case char(9474) % Straight T-B
                temp_top = inOut_T{i-1,j}; % Format input data
                temp_bottom = inOut_T{i+1,j};
                temp_top_C = inOut_C{i-1,j};
                temp_bottom_C = inOut_C{i+1,j};
                temp_top_F = inOut_F{i-1,j};
                temp_bottom_F = inOut_F{i+1,j};
                mlt = [temp_top(1,3);temp(:,1);temp_bottom(2,1)];
                mrt = [temp_bottom(1,1);temp(:,2);temp_top(2,3)];
                mlt_C = [temp_top_C(1,3);temp_C(:,1);temp_bottom_C(2,1)];
                mrt_C = [temp_bottom_C(1,1);temp_C(:,2);temp_top_C(2,3)];
                mlt_F = [temp_top_F(1,3);temp_F(:,1);temp_bottom_F(2,1)];
                mrt_F = [temp_bottom_F(1,1);temp_F(:,2);temp_top_F(2,3)];
                [mlt1,mrt1] = GUI_Straight_CF(mlt,mrt,mlt_C,mrt_C,mlt_F,mrt_F,length);% Calculation function
                inOut_T1(i,j) = {[mrt1(end),0,mlt1(end),0;... % Format output data
                    mlt1(1),0,mrt1(1),0]};
                segments(i,j) = {[mlt1,mrt1]};
            case char(9488) % Corner L-B
                temp_left = inOut_T{i,j-1}; % Format input data
                temp_bottom = inOut_T{i+1,j};
                temp_left_C = inOut_C{i,j-1};
                temp_bottom_C = inOut_C{i+1,j};
                temp_left_F = inOut_F{i,j-1};
                temp_bottom_F = inOut_F{i+1,j};
                mlt = [temp_bottom(1,1);temp(:,1);temp_left(2,2)];
                mrt = [temp_left(1,2);temp(:,2);temp_bottom(2,1)];
                mlt_C = [temp_bottom_C(1,1);temp_C(:,1);temp_left_C(2,2)];
                mrt_C = [temp_left_C(1,2);temp_C(:,2);temp_bottom_C(2,1)];
                mlt_F = [temp_bottom_F(1,1);temp_F(:,1);temp_left_F(2,2)];
                mrt_F = [temp_left_F(1,2);temp_F(:,2);temp_bottom_F(2,1)];
                [mlt1,mrt1] = GUI_Straight_CF(mlt,mrt,mlt_C,mrt_C,mlt_F,mrt_F,length-2); % Calculation function
                inOut_T1(i,j) = {[0,0,mrt1(end),mlt1(end);... % Format output data
                    0,0,mlt1(1),mrt1(1)]};
                segments(i,j) = {[mlt1,mrt1]};
            case char(9484) % Corner B-R
                temp_bottom = inOut_T{i+1,j}; % Format input data
                temp_right = inOut_T{i,j+1};
                temp_bottom_C = inOut_C{i+1,j};
                temp_right_C = inOut_C{i,j+1};
                temp_bottom_F = inOut_F{i+1,j};
                temp_right_F = inOut_F{i,j+1};
                mlt = [temp_right(1,4);temp(:,1);temp_bottom(2,1)];
                mrt = [temp_bottom(1,1);temp(:,2);temp_right(2,4)];
                mlt_C = [temp_right_C(1,4);temp_C(:,1);temp_bottom_C(2,1)];
                mrt_C = [temp_bottom_C(1,1);temp_C(:,2);temp_right_C(2,4)];
                mlt_F = [temp_right_F(1,4);temp_F(:,1);temp_bottom_F(2,1)];
                mrt_F = [temp_bottom_F(1,1);temp_F(:,2);temp_right_F(2,4)];
                [mlt1,mrt1] = GUI_Straight_CF(mlt,mrt,mlt_C,mrt_C,mlt_F,mrt_F,length-2); % Calculation function
                inOut_T1(i,j) = {[0,mrt1(end),mlt1(end),0;... % Format output data
                    0,mlt1(1),mrt1(1),0]};
                segments(i,j) = {[mlt1,mrt1]};
            case char(9492) % Corner R-T
                temp_right = inOut_T{i,j+1}; % Format input data
                temp_top = inOut_T{i-1,j};
                temp_right_C = inOut_C{i,j+1};
                temp_top_C = inOut_C{i-1,j};
                temp_right_F = inOut_F{i,j+1};
                temp_top_F = inOut_F{i-1,j};
                mlt = [temp_top(1,3);temp(:,1);temp_right(2,4)];
                mrt = [temp_right(1,4);temp(:,2);temp_top(2,3)];
                mlt_C = [temp_top_C(1,3);temp_C(:,1);temp_right_C(2,4)];
                mrt_C = [temp_right_C(1,4);temp_C(:,2);temp_top_C(2,3)];
                mlt_F = [temp_top_F(1,3);temp_F(:,1);temp_right_F(2,4)];
                mrt_F = [temp_right_F(1,4);temp_F(:,2);temp_top_F(2,3)];
                [mlt1,mrt1] = GUI_Straight_CF(mlt,mrt,mlt_C,mrt_C,mlt_F,mrt_F,length-2); % Calculation function
                inOut_T1(i,j) = {[mrt1(end),mlt1(end),0,0;... % Format output data
                    mlt1(1),mrt1(1),0,0]};
                segments(i,j) = {[mlt1,mrt1]};
            case char(9496) % Corner T-L
                temp_top = inOut_T{i-1,j}; % Format input data
                temp_left = inOut_T{i,j-1};
                temp_top_C = inOut_C{i-1,j};
                temp_left_C = inOut_C{i,j-1};
                temp_top_F = inOut_F{i-1,j};
                temp_left_F = inOut_F{i,j-1};
                mlt = [temp_left(1,2);temp(:,1);temp_top(2,3)];
                mrt = [temp_top(1,3);temp(:,2);temp_left(2,2)];
                mlt_C = [temp_left_C(1,2);temp_C(:,1);temp_top_C(2,3)];
                mrt_C = [temp_top_C(1,3);temp_C(:,2);temp_left_C(2,2)];
                mlt_F = [temp_left_F(1,2);temp_F(:,1);temp_top_F(2,3)];
                mrt_F = [temp_top_F(1,3);temp_F(:,2);temp_left_F(2,2)];
                [mlt1,mrt1] = GUI_Straight_CF(mlt,mrt,mlt_C,mrt_C,mlt_F,mrt_F,length-2); % Calculation function
                inOut_T1(i,j) = {[mlt1(end),0,0,mrt1(end);... % Format output data
                    mrt1(1),0,0,mlt1(1)]};
                segments(i,j) = {[mlt1,mrt1]};
            case char(9508) % Intersection from Left
                temp_top = inOut_T{i-1,j}; % Format input data
                temp_left = inOut_T{i,j-1};
                temp_bottom = inOut_T{i+1,j};
                temp_top_C = inOut_C{i-1,j};
                temp_left_C = inOut_C{i,j-1};
                temp_bottom_C = inOut_C{i+1,j};
                temp_top_F = inOut_F{i-1,j};
                temp_left_F = inOut_F{i,j-1};
                temp_bottom_F = inOut_F{i+1,j};
                mlt = [temp_top(1,3);temp(1:length,1);temp_bottom(2,1)];
                mrt = [temp_bottom(1,1);temp(1:length,2);temp_top(2,3)];
                slt = [temp_left(1,2);temp(length+1:end-1,1)];
                srt = [temp(length+1:end-1,2);temp_left(2,2)];
                t2t = temp(end,1);
                t1t = temp(end,2);
                mlt_C = [temp_top_C(1,3);temp_C(1:length,1);temp_bottom_C(2,1)];
                mrt_C = [temp_bottom_C(1,1);temp_C(1:length,2);temp_top_C(2,3)];
                slt_C = [temp_left_C(1,2);temp_C(length+1:end-1,1)];
                srt_C = [temp_C(length+1:end-1,2);temp_left_C(2,2)];
                t2t_C = temp_C(end,1);
                t1t_C = temp_C(end,2);
                mlt_F = [temp_top_F(1,3);temp_F(1:length,1);temp_bottom_F(2,1)];
                mrt_F = [temp_bottom_F(1,1);temp_F(1:length,2);temp_top_F(2,3)];
                slt_F = [temp_left_F(1,2);temp_F(length+1:end-1,1)];
                srt_F = [temp_F(length+1:end-1,2);temp_left_F(2,2)];
                t2t_F = temp_F(end,1);
                t1t_F = temp_F(end,2);
                [mrt1,mlt1,srt1,slt1,t1t1,t2t1] = GUI_Intersection_CF(mrt,mlt,srt,slt,t1t,t2t,...
                    mrt_C,mlt_C,srt_C,slt_C,t1t_C,t2t_C,...
                    mrt_F,mlt_F,srt_F,slt_F,t1t_F,t2t_F,length); % Calculation function
                inOut_T1(i,j) = {[mrt1(end),0,mlt1(end),srt1(end);... % Format output data
                    mlt1(1),0,mrt1(1),slt1(1)]};
                segments(i,j) = {[mlt1,mrt1;slt1,srt1;t2t1,t1t1]};
            case char(9516) % Intersection from Bottom
                temp_left = inOut_T{i,j-1}; % Format input data
                temp_bottom = inOut_T{i+1,j};
                temp_right = inOut_T{i,j+1};
                temp_left_C = inOut_C{i,j-1};
                temp_bottom_C = inOut_C{i+1,j};
                temp_right_C = inOut_C{i,j+1};
                temp_left_F = inOut_F{i,j-1};
                temp_bottom_F = inOut_F{i+1,j};
                temp_right_F = inOut_F{i,j+1};
                mlt = [temp_left(1,2);temp(1:length,1);temp_right(2,4)];
                mrt = [temp_right(1,4);temp(1:length,2);temp_left(2,2)];
                slt = [temp_bottom(1,1);temp(length+1:end-1,1)];
                srt = [temp(length+1:end-1,2);temp_bottom(2,1)];
                t2t = temp(end,1);
                t1t = temp(end,2);
                mlt_C = [temp_left_C(1,2);temp_C(1:length,1);temp_right_C(2,4)];
                mrt_C = [temp_right_C(1,4);temp_C(1:length,2);temp_left_C(2,2)];
                slt_C = [temp_bottom_C(1,1);temp_C(length+1:end-1,1)];
                srt_C = [temp_C(length+1:end-1,2);temp_bottom_C(2,1)];
                t2t_C = temp_C(end,1);
                t1t_C = temp_C(end,2);
                mlt_F = [temp_left_F(1,2);temp_F(1:length,1);temp_right_F(2,4)];
                mrt_F = [temp_right_F(1,4);temp_F(1:length,2);temp_left_F(2,2)];
                slt_F = [temp_bottom_F(1,1);temp_F(length+1:end-1,1)];
                srt_F = [temp_F(length+1:end-1,2);temp_bottom_F(2,1)];
                t2t_F = temp_F(end,1);
                t1t_F = temp_F(end,2);
                [mrt1,mlt1,srt1,slt1,t1t1,t2t1] = GUI_Intersection_CF(mrt,mlt,srt,slt,t1t,t2t,...
                    mrt_C,mlt_C,srt_C,slt_C,t1t_C,t2t_C,...
                    mrt_F,mlt_F,srt_F,slt_F,t1t_F,t2t_F,length); % Calculation function
                inOut_T1(i,j) = {[0,mlt1(end),srt1(end),mrt1(end);... % Format output data
                    0,mrt1(1),slt1(1),mlt1(1)]};
                segments(i,j) = {[mlt1,mrt1;slt1,srt1;t2t1,t1t1]};
            case char(9500) % Intersection from Right
                temp_top = inOut_T{i-1,j}; % Format input data
                temp_right = inOut_T{i,j+1};
                temp_bottom = inOut_T{i+1,j};
                temp_top_C = inOut_C{i-1,j};
                temp_right_C = inOut_C{i,j+1};
                temp_bottom_C = inOut_C{i+1,j};
                temp_top_F = inOut_F{i-1,j};
                temp_right_F = inOut_F{i,j+1};
                temp_bottom_F = inOut_F{i+1,j};
                mlt = [temp_bottom(1,1);temp(1:length,1);temp_top(2,3)];
                mrt = [temp_top(1,3);temp(1:length,2);temp_bottom(2,1)];
                slt = [temp_right(1,4);temp(length+1:end-1,1)];
                srt = [temp(length+1:end-1,2);temp_right(2,4)];
                t2t = temp(end,1);
                t1t = temp(end,2);
                mlt_C = [temp_bottom_C(1,1);temp_C(1:length,1);temp_top_C(2,3)];
                mrt_C = [temp_top_C(1,3);temp_C(1:length,2);temp_bottom_C(2,1)];
                slt_C = [temp_right_C(1,4);temp_C(length+1:end-1,1)];
                srt_C = [temp_C(length+1:end-1,2);temp_right_C(2,4)];
                t2t_C = temp_C(end,1);
                t1t_C = temp_C(end,2);
                mlt_F = [temp_bottom_F(1,1);temp_F(1:length,1);temp_top_F(2,3)];
                mrt_F = [temp_top_F(1,3);temp_F(1:length,2);temp_bottom_F(2,1)];
                slt_F = [temp_right_F(1,4);temp_F(length+1:end-1,1)];
                srt_F = [temp_F(length+1:end-1,2);temp_right_F(2,4)];
                t2t_F = temp_F(end,1);
                t1t_F = temp_F(end,2);
                [mrt1,mlt1,srt1,slt1,t1t1,t2t1] = GUI_Intersection_CF(mrt,mlt,srt,slt,t1t,t2t,...
                    mrt_C,mlt_C,srt_C,slt_C,t1t_C,t2t_C,...
                    mrt_F,mlt_F,srt_F,slt_F,t1t_F,t2t_F,length); % Calculation function
                inOut_T1(i,j) = {[mlt1(end),srt1(end),mrt1(end),0;... % Format output data
                    mrt1(1),slt1(1),mlt1(1),0]};
                segments(i,j) = {[mlt1,mrt1;slt1,srt1;t2t1,t1t1]};
            case char(9524) % Intersection from Top
                temp_top = inOut_T{i-1,j}; % Format input data
                temp_left = inOut_T{i,j-1};
                temp_right = inOut_T{i,j+1};
                temp_top_C = inOut_C{i-1,j};
                temp_left_C = inOut_C{i,j-1};
                temp_right_C = inOut_C{i,j+1};
                temp_top_F = inOut_F{i-1,j};
                temp_left_F = inOut_F{i,j-1};
                temp_right_F = inOut_F{i,j+1};
                mlt = [temp_right(1,4);temp(1:length,1); temp_left(2,2)];
                mrt = [temp_left(1,2);temp(1:length,2);temp_right(2,4)];
                slt = [temp_top(1,3);temp(length+1:end-1,1)];
                srt = [temp(length+1:end-1,2);temp_top(2,3)];
                t2t = temp(end,1);
                t1t = temp(end,2);
                mlt_C = [temp_right_C(1,4);temp_C(1:length,1);temp_left_C(2,2)];
                mrt_C = [temp_left_C(1,2);temp_C(1:length,2);temp_right_C(2,4)];
                slt_C = [temp_top_C(1,3);temp_C(length+1:end-1,1)];
                srt_C = [temp_C(length+1:end-1,2);temp_top_C(2,3)];
                t2t_C = temp_C(end,1);
                t1t_C = temp_C(end,2);
                mlt_F = [temp_right_F(1,4);temp_F(1:length,1);temp_left_F(2,2)];
                mrt_F = [temp_left_F(1,2);temp_F(1:length,2);temp_right_F(2,4)];
                slt_F = [temp_top_F(1,3);temp_F(length+1:end-1,1)];
                srt_F = [temp_F(length+1:end-1,2);temp_top_F(2,3)];
                t2t_F = temp_F(end,1);
                t1t_F = temp_F(end,2);
                [mrt1,mlt1,srt1,slt1,t1t1,t2t1] = GUI_Intersection_CF(mrt,mlt,srt,slt,t1t,t2t,...
                    mrt_C,mlt_C,srt_C,slt_C,t1t_C,t2t_C,...
                    mrt_F,mlt_F,srt_F,slt_F,t1t_F,t2t_F,length); % Calculation function
                inOut_T1(i,j) = {[srt1(end),mrt1(end),0,mlt1(end);... % Format output data
                    slt1(1),mlt1(1),0,mrt1(1)]};
                segments(i,j) = {[mlt1,mrt1;slt1,srt1;t2t1,t1t1]};
            case char(9711) % Roundabout
                temp_top = inOut_T{i-1,j}; % Format input data
                temp_right = inOut_T{i,j+1};
                temp_bottom = inOut_T{i+1,j};
                temp_left = inOut_T{i,j-1};
                temp_top_C = inOut_C{i-1,j};
                temp_right_C = inOut_C{i,j+1};
                temp_bottom_C = inOut_C{i+1,j};
                temp_left_C = inOut_C{i,j-1};
                temp_top_F = inOut_F{i-1,j};
                temp_right_F = inOut_F{i,j+1};
                temp_bottom_F = inOut_F{i+1,j};
                temp_left_F = inOut_F{i,j-1};
                mtt = [temp(end-2,2);temp(1:end-2,1);temp(1,4)]; % Top main road
                mrt = [temp(end-2,3);temp(1:end-2,2);temp(1,1)]; % Right main road
                mbt = [temp(end-2,4);temp(1:end-2,3);temp(1,2)]; % Bottom main road
                mlt = [temp(end-2,1);temp(1:end-2,4);temp(1,3)]; % Left main road
                stlt = [temp_top(1,3);temp(end-1,1)]; % Left side road of top portion of roundabout
                strt = [temp(end,1);temp_top(2,3)]; % Right side road of top portion of roundabout
                srlt = [temp_right(1,4);temp(end-1,2)];
                srrt = [temp(end,2);temp_right(2,4)];
                sblt = [temp_bottom(1,1);temp(end-1,3)];
                sbrt = [temp(end,3);temp_bottom(2,1)];
                sllt = [temp_left(1,2);temp(end-1,4)];
                slrt = [temp(end,4);temp_right(2,2)];
                mtt_C = [temp_C(end-2,2);temp_C(1:end-2,1);temp_C(1,4)]; % Top main road
                mrt_C = [temp_C(end-2,3);temp_C(1:end-2,2);temp_C(1,1)]; % Right main road
                mbt_C = [temp_C(end-2,4);temp_C(1:end-2,3);temp_C(1,2)]; % Bottom main road
                mlt_C = [temp_C(end-2,1);temp_C(1:end-2,4);temp_C(1,3)]; % Left main road
                stlt_C = [temp_top_C(1,3);temp_C(end-1,1)]; % Left side road of top portion of roundabout
                strt_C = [temp_C(end,1);temp_top_C(2,3)]; % Right side road of top portion of roundabout
                srlt_C = [temp_right_C(1,4);temp_C(end-1,2)];
                srrt_C = [temp_C(end,2);temp_right_C(2,4)];
                sblt_C = [temp_bottom_C(1,1);temp_C(end-1,3)];
                sbrt_C = [temp_C(end,3);temp_bottom_C(2,1)];
                sllt_C = [temp_left_C(1,2);temp_C(end-1,4)];
                slrt_C = [temp_C(end,4);temp_right_C(2,2)];
                mtt_F = [temp_F(end-2,2);temp_F(1:end-2,1);temp_F(1,4)]; % Top main road
                mrt_F = [temp_F(end-2,3);temp_F(1:end-2,2);temp_F(1,1)]; % Right main road
                mbt_F = [temp_F(end-2,4);temp_F(1:end-2,3);temp_F(1,2)]; % Bottom main road
                mlt_F = [temp_F(end-2,1);temp_F(1:end-2,4);temp_F(1,3)]; % Left main road
                stlt_F = [temp_top_F(1,3);temp_F(end-1,1)]; % Left side road of top portion of roundabout
                strt_F = [temp_F(end,1);temp_top_F(2,3)]; % Right side road of top portion of roundabout
                srlt_F = [temp_right_F(1,4);temp_F(end-1,2)];
                srrt_F = [temp_F(end,2);temp_right_F(2,4)];
                sblt_F = [temp_bottom_F(1,1);temp_F(end-1,3)];
                sbrt_F = [temp_F(end,3);temp_bottom_F(2,1)];
                sllt_F = [temp_left_F(1,2);temp_F(end-1,4)];
                slrt_F = [temp_F(end,4);temp_right_F(2,2)];
                [mtt1,mrt1,mbt1,mlt1,stlt1,strt1,srlt1,srrt1,sblt1,sbrt1,sllt1,slrt1]...
                    = GUI_Roundabout_CF(mtt,mrt,mbt,mlt,stlt,strt,srlt,srrt,sblt,sbrt,sllt,slrt,...
                    mtt_C,mrt_C,mbt_C,mlt_C,stlt_C,strt_C,srlt_C,srrt_C,sblt_C,sbrt_C,sllt_C,slrt_C,...
                    mtt_F,mrt_F,mbt_F,mlt_F,stlt_F,strt_F,srlt_F,srrt_F,sblt_F,sbrt_F,sllt_F,slrt_F,length); % Calculation function
                inOut_T1(i,j) = {[strt1,srrt1,sbrt1,slrt1;... % Format output data
                    stlt1,srlt1,sblt1,sllt1]};
                segments(i,j) = {[mtt1,mrt1,mbt1,mlt1;stlt1,srlt1,sblt1,sllt1;strt1,srrt1,sbrt1,slrt1]};
            case char(9532) % 4-Way traffic light
                ls = length/2-3/2; % Length of the side roads
                temp_top = inOut_T{i-1,j}; % Format input data
                temp_right = inOut_T{i,j+1};
                temp_bottom = inOut_T{i+1,j};
                temp_left = inOut_T{i,j-1};
                temp_top_C = inOut_C{i-1,j};
                temp_right_C = inOut_C{i,j+1};
                temp_bottom_C = inOut_C{i+1,j};
                temp_left_C = inOut_C{i,j-1};
                temp_top_F = inOut_F{i-1,j};
                temp_right_F = inOut_F{i,j+1};
                temp_bottom_F = inOut_F{i+1,j};
                temp_left_F = inOut_F{i,j-1};
                tl = [temp_top(1,3);temp(1:ls,1)];
                tr = [temp(ls+1:end-2,1);temp_top(2,3)];
                rl = [temp_right(1,4);temp(1:ls,2)];
                rr = [temp(ls+1:end-2,2);temp_right(2,4)];
                bl = [temp_bottom(1,1);temp(1:ls,3)];
                br = [temp(ls+1:end-2,3);temp_bottom(2,1)];
                ll = [temp_left(1,2);temp(1:ls,4)];
                lr = [temp(ls+1:end-2,4);temp_left(2,2)];
                tc = temp(end-1,1);
                rc = temp(end-1,2);
                bc = temp(end-1,3);
                lc = temp(end-1,4);
                tl_C = [temp_top_C(1,3);temp_C(1:ls,1)];
                tr_C = [temp_C(ls+1:end-1,1);temp_top_C(2,3)];
                rl_C = [temp_right_C(1,4);temp_C(1:ls,2)];
                rr_C = [temp_C(ls+1:end-1,2);temp_right_C(2,4)];
                bl_C = [temp_bottom_C(1,1);temp_C(1:ls,3)];
                br_C = [temp_C(ls+1:end-1,3);temp_bottom_C(2,1)];
                ll_C = [temp_left_C(1,2);temp_C(1:ls,4)];
                lr_C = [temp_C(ls+1:end-1,4);temp_left_C(2,2)];
                tc_C = temp_C(end,1);
                rc_C = temp_C(end,2);
                bc_C = temp_C(end,3);
                lc_C = temp_C(end,4);
                tl_F = [temp_top_F(1,3);temp_F(1:ls,1)];
                tr_F = [temp_F(ls+1:end-1,1);temp_top_F(2,3)];
                rl_F = [temp_right_F(1,4);temp_F(1:ls,2)];
                rr_F = [temp_F(ls+1:end-1,2);temp_right_F(2,4)];
                bl_F = [temp_bottom_F(1,1);temp_F(1:ls,3)];
                br_F = [temp_F(ls+1:end-1,3);temp_bottom_F(2,1)];
                ll_F = [temp_left_F(1,2);temp_F(1:ls,4)];
                lr_F = [temp_F(ls+1:end-1,4);temp_left_F(2,2)];
                tc_F = temp_F(end,1);
                rc_F = temp_F(end,2);
                bc_F = temp_F(end,3);
                lc_F = temp_F(end,4);
                timer = temp(end,1);
                [tl_1,tr_1,rl_1,rr_1,bl_1,br_1,ll_1,lr_1,tc_1,rc_1,bc_1,lc_1,timer]...
                    = GUI_Traffic_Light_CF(tl,tr,rl,rr,bl,br,ll,lr,tc,rc,bc,lc,...
                    tl_C,tr_C,rl_C,rr_C,bl_C,br_C,ll_C,lr_C,tc_C,rc_C,bc_C,lc_C,...
                    tl_F,tr_F,rl_F,rr_F,bl_F,br_F,ll_F,lr_F,tc_F,rc_F,bc_F,lc_F,length,timer); % Calculation function
                inOut_T1(i,j) = {[tr_1(end),rr_1(end),br_1(end),lr_1(end);... % Format output data
                    tl_1(1),rl_1(1),bl_1(1),ll_1(1)]};
                segments(i,j) = {[tl_1,rl_1,bl_1,ll_1;tr_1,rr_1,br_1,lr_1;tc_1,rc_1,bc_1,lc_1;timer,0,0,0]};
            otherwise % Case of '', not used cell        
                
        end
        
    end
end

handles.Segments = segments;
inOut_T = inOut_T1;
handles.InOut_T = inOut_T;
handles.InOut_T1 = inOut_T1;

visualize_Network(handles); % Visualize results

guidata(hObject,handles); % Update global variables

function visualize_Network(handles)

% Retrieve global variables
length = handles.SegmentLength;
simSize = handles.Size;
road_Network = handles.Network;
segments = handles.Segments;
inOut_T = handles.InOut_T;

visualization = zeros(length*simSize);
vis_Segment = zeros(length); % Segment used to format the stored data for visualization

cl = length/2-0.5; % Constant, position of left lane
cr = length/2+1.5; % Constant, position of right lane

for i = 1:simSize
    for j = 1:simSize
        
        vis_Segment(:) = 0;
        
        temp = segments{i,j}+0.05;
        
        switch road_Network(i,j)
            case char(9472) % Straight L-R
                vis_Segment(cl,:) = flipud(temp(:,1));
                vis_Segment(cr,:) = temp(:,2);
            case char(9474) % Straight T-B
                vis_Segment(:,cl) = temp(:,1);
                vis_Segment(:,cr) = flipud(temp(:,2));
            case char(9488) % Corner L-B
                vis_Segment(cr,1:cl) = temp(1:cl,2);
                vis_Segment(cr+1:end,cl) = temp(cl+1:end,2);
                vis_Segment(cr+2:end,cr) = flipud(temp(1:cl-2,1));
                vis_Segment(cr:cr+1,cr) = flipud(temp(cl-1,1));
                vis_Segment(cl:cl+1,cr) = flipud(temp(cl,1));
                vis_Segment(cl,cl+1) = flipud(temp(cl,1));
                vis_Segment(cl,cl-1:cl) = flipud(temp(cl+1,1));
                vis_Segment(cl,1:cl-2) = flipud(temp(cl+2:end,1));
            case char(9484) % Corner B-R
                vis_Segment(cr:end,cr) = flipud(temp(1:cl,2));
                vis_Segment(cr,cr+1:end) = temp(cl+1:end,2);
                vis_Segment(cr+2:end,cl) = temp(cl+2:end,1);
                vis_Segment(cr:cr+1,cl) = temp(cl+1,1);
                vis_Segment(cl:cl+1,cl) = temp(cl,1);
                vis_Segment(cl,cl+1) = temp(cl,1);
                vis_Segment(cl,cr:cr+1) = temp(cl-1,1);
                vis_Segment(cl,cr+2:end) = flipud(temp(1:cl-2,1));
            case char(9492) % Corner R-T
                vis_Segment(cl,cr:end) = flipud(temp(1:cl,2));
                vis_Segment(1:cl-1,cr) = flipud(temp(cl+1:end,2));
                vis_Segment(1:cl-2,cl) = temp(1:cl-2,1);
                vis_Segment(cl-1:cl,cl) = temp(cl-1,1);
                vis_Segment(cl+1:cl+2,cl) = temp(cl,1);
                vis_Segment(cr,cl+1) = temp(cl,1);
                vis_Segment(cr,cr:cr+1) = temp(cl+1,1);
                vis_Segment(cr,cr+2:end) = temp(cl+2:end,1);
            case char(9496) % Corner T-L
                vis_Segment(1:cl,cl) = temp(1:cl,2);
                vis_Segment(cl,1:cl-1) = flipud(temp(cl+1:end,2));
                vis_Segment(cr,1:cl-2) = temp(1:cl-2,1);
                vis_Segment(cr,cl-1:cl) = temp(cl-1,1);
                vis_Segment(cr,cl+1:cl+2) = temp(cl,1);
                vis_Segment(cl+1,cr) = temp(cl,1);
                vis_Segment(cl-1:cl,cr) = temp(cl+1,1);
                vis_Segment(1:cl-2,cr) = flipud(temp(cl+2:end,1));
            case char(9508) % Intersection Left
                vis_Segment(:,cl) = temp(1:length,1);
                vis_Segment(:,cr) = flipud(temp(1:length,2));
                vis_Segment(cl,1:cl-1) = flipud(temp(length+1:end-1,2));
                vis_Segment(cr,1:cl-1) = temp(length+1:end-1,1);
                vis_Segment(cl,cl+1) = temp(end,2);
                vis_Segment(cr,cl+1) = temp(end,1);
            case char(9516) % Intersection Bottom
                vis_Segment(cl,:) = flipud(temp(1:length,2));
                vis_Segment(cr,:) = temp(1:length,1);
                vis_Segment(cr+1:end,cl) = temp(length+1:end-1,2);
                vis_Segment(cr+1:end,cr) = flipud(temp(length+1:end-1,1));
                vis_Segment(cl+1,cl) = temp(end,2);
                vis_Segment(cl+1,cr) = temp(end,1);
            case char(9500) % Intersection Right
                vis_Segment(:,cl) = temp(1:length,2);
                vis_Segment(:,cr) = flipud(temp(1:length,1));
                vis_Segment(cr,cr+1:end) = temp(length+1:end-1,2);
                vis_Segment(cl,cr+1:end) = flipud(temp(length+1:end-1,1));
                vis_Segment(cr,cl+1) = temp(end,2);
                vis_Segment(cl,cl+1) = temp(end,1);
            case char(9524) % Intersection Top
                %temp = temp-0.2;
                vis_Segment(cr,:) = temp(1:length,2);
                vis_Segment(cl,:) = flipud(temp(1:length,1));
                vis_Segment(1:cl-1,cr) = flipud(temp(length+1:end-1,2));
                vis_Segment(1:cl-1,cl) = temp(length+1:end-1,1);
                vis_Segment(cl+1,cr) = temp(end,2);
                vis_Segment(cl+1,cl) = temp(end,1);
                %disp(vis_Segment);
            case char(9711) % Roundabout
                vis_Segment(2,3:end-1) = flipud(temp(1:end-2,1));
                vis_Segment(3:end-1,end-1) = flipud(temp(1:end-2,2));
                vis_Segment(end-1,2:end-2) = temp(1:end-2,3);
                vis_Segment(2:end-2,2) = temp(1:end-2,4);
                vis_Segment(1,cl) = temp(end-1,1);
                vis_Segment(1,cr) = temp(end,1);
                vis_Segment(cl,end) = temp(end-1,2);
                vis_Segment(cr,end) = temp(end,2);
                vis_Segment(end,cr) = temp(end-1,3);
                vis_Segment(end,cl) = temp(end,3);
                vis_Segment(cr,1) = temp(end-1,4);
                vis_Segment(cl,1) = temp(end,4);
            case char(9532)
                ls = length/2-3/2; % Length of the side roads
                tl = temp(1:ls,1);
                tr = temp(ls+1:end-2,1);
                rl = temp(1:ls,2);
                rr = temp(ls+1:end-2,2);
                bl = temp(1:ls,3);
                br = temp(ls+1:end-2,3);
                ll = temp(1:ls,4);
                lr = temp(ls+1:end-2,4);
                tc = temp(end-1,1);
                rc = temp(end-1,2);
                bc = temp(end-1,3);
                lc = temp(end-1,4);
                vis_Segment(1:ls,cl) = tl;
                vis_Segment(cl,end-(ls-1):end) = flipud(rl);
                vis_Segment(end-(ls-1):end,cr) = flipud(bl);
                vis_Segment(cr,1:ls) = ll;
                vis_Segment(1:ls+1,cr) = flipud(tr);
                vis_Segment(cr,end-ls:end) = rr;
                vis_Segment(end-ls:end,cl) = br;
                vis_Segment(cl,1:ls+1) = flipud(lr);
                vis_Segment(cl,cl+1) = tc;
                vis_Segment(cl+1,cr) = rc;
                vis_Segment(cr,cl+1) = bc;
                vis_Segment(cl+1,cl) = lc;
            case 'EX'
                temp = inOut_T{i,j}+0.2;
                if i-1 == 0 % Top side
                elseif ~strcmp(road_Network(i-1,j),'')
                    vis_Segment(1,cr) = temp(1,1);
                end
                if j+1 > simSize % Right side
                elseif ~strcmp(road_Network(i,j+1),'')
                    vis_Segment(cr,end) = temp(1,2);
                end
                if i+1 > simSize % Bottom side
                elseif ~strcmp(road_Network(i+1,j),'')
                    vis_Segment(end,cl) = temp(1,3);
                end
                if j-1 == 0 % Left side
                elseif ~strcmp(road_Network(i,j-1),'')
                    vis_Segment(cl,1) = temp(1,4);
                end
            otherwise % Case of empty cell ''
                
        end
        
        visualization(1+(i-1)*length:i*length,1+(j-1)*length:j*length) = vis_Segment; % Add visualization segment to overall visualization matrix
        
    end
end
heatmap(visualization,'Colormap',handles.cMap,'ColorLimits',[0,3.05],'GridVisible','off','CellLabelColor','none'); %Create heatmap of data

drawnow;



function save_Network(handles)

%Retrieve global variables
fileName = handles.FileName;
length = handles.SegmentLength;
simSize = handles.Size;
road_Network = handles.Network;
segments = handles.Segments;
inOut_T = handles.InOut_T;
inOut_T1 = handles.InOut_T1;
inputVolume = handles.InputVolume;
connections = handles.Connections;

fileID = fopen(fileName,'w'); % Open file
fprintf(fileID,'%i\n%i\n',simSize,length); % Save size of system, and length of individual segments
if strcmp(inputVolume,'random')
    fprintf(fileID,[inputVolume newline]); % Save the input volume, rate at which vehicles enter the system
else
    fprintf(fileID,'%i\n',inputVolume);
end

for i = 1:simSize % Loop over all segments
    
    for j = 1:simSize  
        
        temp_Segment = segments{i,j}; % Retrieve data for the segment currently being saved
        temp_T = inOut_T{i,j}; 
        temp_T1 = inOut_T1{i,j};
        temp_Connections = connections{i,j};
        fprintf(fileID,'%i',char(road_Network(i,j))); % Save the type of road segment
        fprintf(fileID,newline);
        fprintf(fileID,'%d ',temp_Segment(:)); % Save cell data
        fprintf(fileID, newline);
        fprintf(fileID,'%d ',temp_T(:)); % In- and output data for this and last time step
        fprintf(fileID, newline);
        fprintf(fileID,'%d ',temp_T1(:));
        fprintf(fileID, newline);
        fprintf(fileID,'%i ',temp_Connections(:)); % Save connections matrix of the segment
        fprintf(fileID, newline);
        
    end
    
end

fclose(fileID); % Close file



function load_Network(hObject,handles)

% Load global variables
fileName = handles.FileName;

fileID = fopen(fileName,'r');

handles.Size = str2double(fgetl(fileID)); % Retrieve size of system
segLength = str2double(fgetl(fileID)); % Retrieve length of segments
handles.SegmentLength = segLength;
guidata(hObject,handles); % Save global variables

road_Network = strings(handles.Size); % Initialize the global variables for system being loaded
connections = cell(handles.Size);
segments = cell(handles.Size);
inOut_T = cell(handles.Size);
inOut_T1 = cell(handles.Size);

temp = fgetl(fileID);
if strcmp(temp,'random') % Retrieve input volume
    handles.InputVolume = 'random';
else
    handles.InputVolume = str2double(temp);
end

for i = 1:handles.Size % Loop over all imported segments
    
    for j = 1:handles.Size
        
        temp_Type = fgetl(fileID);
        
        switch temp_Type
            
            case '6988' % Entrance and exit
                road_Network(i,j) = 'EX'; % Retrieve segment type
                fgetl(fileID);
                segments(i,j) = {0}; % Retrieve segment data
                inOut_T(i,j) = {reshape(str2num(fgetl(fileID)),2,4)}; % Retrieve in- and output data for segment
                inOut_T1(i,j) = {reshape(str2num(fgetl(fileID)),2,4)};
                connections(i,j) = {reshape(str2num(fgetl(fileID)),1,4)}; % Retrieve connections data for segment
            case '' % Empty cell
                road_Network(i,j) = '';
                fgetl(fileID);
                segments(i,j) = {0};
                fgetl(fileID);
                inOut_T(i,j) = {[]}; 
                fgetl(fileID);
                inOut_T1(i,j) = {[]};
                connections(i,j) = {reshape(str2num(fgetl(fileID)),1,4)};
            case {'9472', '9474'} % Straight road segment
                road_Network(i,j) = char(str2num(temp_Type));
                segments(i,j) = {reshape(str2num(fgetl(fileID)),segLength,2)};
                inOut_T(i,j) = {reshape(str2num(fgetl(fileID)),2,4)};
                inOut_T1(i,j) = {reshape(str2num(fgetl(fileID)),2,4)};
                connections(i,j) = {reshape(str2num(fgetl(fileID)),1,4)};
            case {'9488', '9484', '9492', '9496'} % Corner road segment
                road_Network(i,j) = char(str2num(temp_Type));
                segments(i,j) = {reshape(str2num(fgetl(fileID)),segLength-2,2)};
                inOut_T(i,j) = {reshape(str2num(fgetl(fileID)),2,4)};
                inOut_T1(i,j) = {reshape(str2num(fgetl(fileID)),2,4)};
                connections(i,j) = {reshape(str2num(fgetl(fileID)),1,4)};
            case {'9508', '9516', '9500', '9524'} % Intersections
                road_Network(i,j) = char(str2num(temp_Type));
                segments(i,j) = {reshape(str2num(fgetl(fileID)),segLength+segLength/2-1/2,2)};
                inOut_T(i,j) = {reshape(str2num(fgetl(fileID)),2,4)};
                inOut_T1(i,j) = {reshape(str2num(fgetl(fileID)),2,4)};
                connections(i,j) = {reshape(str2num(fgetl(fileID)),1,4)};
            case '9711' % Roundabouts
                road_Network(i,j) = char(str2num(temp_Type));
                segments(i,j) = {reshape(str2num(fgetl(fileID)),segLength-1,4)};
                inOut_T(i,j) = {reshape(str2num(fgetl(fileID)),2,4)};
                inOut_T1(i,j) = {reshape(str2num(fgetl(fileID)),2,4)};
                connections(i,j) = {reshape(str2num(fgetl(fileID)),1,4)};
            case '9532' % Traffic Light
                road_Network(i,j) = char(str2num(temp_Type));
                segments(i,j) = {reshape(str2num(fgetl(fileID)),segLength,4)};
                inOut_T(i,j) = {reshape(str2num(fgetl(fileID)),2,4)};
                inOut_T1(i,j) = {reshape(str2num(fgetl(fileID)),2,4)};
                connections(i,j) = {reshape(str2num(fgetl(fileID)),1,4)};
        end
        
    end
    
end

fclose(fileID);

handles.Network = road_Network;
handles.Segments = segments;
handles.InOut_T = inOut_T;
handles.InOut_T1 = inOut_T1;
handles.Connections = connections;

guidata(hObject,handles); % Save global variables

visualize_Network(handles); % Visualize loaded system
