function varargout = CapandFlow_Edit(varargin)
% CAPANDFLOW_EDIT MATLAB code for CapandFlow_Edit.fig
%      CAPANDFLOW_EDIT, by itself, creates a new CAPANDFLOW_EDIT or raises the existing
%      singleton*.
%
%      H = CAPANDFLOW_EDIT returns the handle to a new CAPANDFLOW_EDIT or the handle to
%      the existing singleton*.
%
%      CAPANDFLOW_EDIT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CAPANDFLOW_EDIT.M with the given input arguments.
%
%      CAPANDFLOW_EDIT('Property','Value',...) creates a new CAPANDFLOW_EDIT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before CapandFlow_Edit_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to CapandFlow_Edit_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help CapandFlow_Edit

% Last Modified by GUIDE v2.5 13-Jul-2019 00:20:48

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @CapandFlow_Edit_OpeningFcn, ...
                   'gui_OutputFcn',  @CapandFlow_Edit_OutputFcn, ...
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


% --- Executes just before CapandFlow_Edit is made visible.
function CapandFlow_Edit_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to CapandFlow_Edit (see VARARGIN)

% Choose default command line output for CapandFlow_Edit
handles.output = hObject;

if nargin<4 | ~isstruct(varargin{1})
else
    
    road_Network = varargin{1}.Network;
    CorF = varargin{1}.CorF;
    Segi = varargin{1}.Segi;
    Segj = varargin{1}.Segj;
    connections = varargin{2};
    segments = varargin{3};
    segLength = varargin{1}.segmentLength;
    capacity = varargin{4};
    flow = varargin{5};
    CPWindow = varargin{6};
    dist = 0.7/(segLength); % Variables for the arrangment of the boxes
    boxSize = 0.7/(segLength);

    handles.Network = road_Network; 
    handles.Connections = connections;
    handles.Segments = segments;
    handles.segmentLength = segLength;
    handles.CPWindow = CPWindow;
    handles.Capacity = capacity;
    handles.Flow = flow;
    handles.CorF = CorF;
    handles.Segi = Segi;
    handles.Segj = Segj;
    
    if strcmp(CorF,'Capacity')
        data = capacity{Segi,Segj};
    else
        data = flow{Segi,Segj};
    end


    vis_Segment = zeros(segLength); % Segment used to format the stored data for visualization

    cl = segLength/2-0.5; % Constant, position of left lane
    cr = segLength/2+1.5; % Constant, position of right lane
        
    vis_Segment(:) =NaN;

    switch road_Network(Segi,Segj)
        case char(9472) % Straight L-R
            vis_Segment(cl,:) = flipud(data(:,1));
            vis_Segment(cr,:) = data(:,2);
        case char(9474) % Straight T-B
            vis_Segment(:,cl) = data(:,1);
            vis_Segment(:,cr) = flipud(data(:,2));
        case char(9488) % Corner L-B
            vis_Segment(cr,1:cl) = data(1:cl,2);
            vis_Segment(cr+1:end,cl) = data(cl+1:end,2);
            vis_Segment(cr+2:end,cr) = flipud(data(1:cl-2,1));
            vis_Segment(cr:cr+1,cr) = flipud(data(cl-1,1));
            vis_Segment(cl:cl+1,cr) = flipud(data(cl,1));
            vis_Segment(cl,cl+1) = flipud(data(cl,1));
            vis_Segment(cl,cl-1:cl) = flipud(data(cl+1,1));
            vis_Segment(cl,1:cl-2) = flipud(data(cl+2:end,1));
        case char(9484) % Corner B-R
            vis_Segment(cr:end,cr) = flipud(data(1:cl,2));
            vis_Segment(cr,cr+1:end) = data(cl+1:end,2);
            vis_Segment(cr+2:end,cl) = data(cl+2:end,1);
            vis_Segment(cr:cr+1,cl) = data(cl+1,1);
            vis_Segment(cl:cl+1,cl) = data(cl,1);
            vis_Segment(cl,cl+1) = data(cl,1);
            vis_Segment(cl,cr:cr+1) = data(cl-1,1);
            vis_Segment(cl,cr+2:end) = flipud(data(1:cl-2,1));
        case char(9492) % Corner R-T
            vis_Segment(cl,cr:end) = flipud(data(1:cl,2));
            vis_Segment(1:cl-1,cr) = flipud(data(cl+1:end,2));
            vis_Segment(1:cl-2,cl) = data(1:cl-2,1);
            vis_Segment(cl-1:cl,cl) = data(cl-1,1);
            vis_Segment(cl+1:cl+2,cl) = data(cl,1);
            vis_Segment(cr,cl+1) = data(cl,1);
            vis_Segment(cr,cr:cr+1) = data(cl+1,1);
            vis_Segment(cr,cr+2:end) = data(cl+2:end,1);
        case char(9496) % Corner T-L
            vis_Segment(1:cl,cl) = data(1:cl,2);
            vis_Segment(cl,1:cl-1) = flipud(data(cl+1:end,2));
            vis_Segment(cr,1:cl-2) = data(1:cl-2,1);
            vis_Segment(cr,cl-1:cl) = data(cl-1,1);
            vis_Segment(cr,cl+1:cl+2) = data(cl,1);
            vis_Segment(cl+1,cr) = data(cl,1);
            vis_Segment(cl-1:cl,cr) = data(cl+1,1);
            vis_Segment(1:cl-2,cr) = flipud(data(cl+2:end,1));
        case char(9508) % Intersection Left
            vis_Segment(:,cl) = data(1:segLength,1);
            vis_Segment(:,cr) = flipud(data(1:segLength,2));
            vis_Segment(cl,1:cl-1) = flipud(data(segLength+1:end-1,2));
            vis_Segment(cr,1:cl-1) = data(segLength+1:end-1,1);
            vis_Segment(cl,cl+1) = data(end,2);
            vis_Segment(cr,cl+1) = data(end,1);
        case char(9516) % Intersection Bottom
            vis_Segment(cl,:) = flipud(data(1:segLength,2));
            vis_Segment(cr,:) = data(1:segLength,1);
            vis_Segment(cr+1:end,cl) = data(segLength+1:end-1,2);
            vis_Segment(cr+1:end,cr) = flipud(data(segLength+1:end-1,1));
            vis_Segment(cl+1,cl) = data(end,2);
            vis_Segment(cl+1,cr) = data(end,1);
        case char(9500) % Intersection Right
            vis_Segment(:,cl) = data(1:segLength,2);
            vis_Segment(:,cr) = flipud(data(1:segLength,1));
            vis_Segment(cr,cr+1:end) = data(segLength+1:end-1,2);
            vis_Segment(cl,cr+1:end) = flipud(data(segLength+1:end-1,1));
            vis_Segment(cr,cl+1) = data(end,2);
            vis_Segment(cl,cl+1) = data(end,1);
        case char(9524) % Intersection Top
            vis_Segment(cr,:) = data(1:segLength,2);
            vis_Segment(cl,:) = flipud(data(1:segLength,1));
            vis_Segment(1:cl-1,cr) = flipud(data(segLength+1:end-1,2));
            vis_Segment(1:cl-1,cl) = data(segLength+1:end-1,1);
            vis_Segment(cl+1,cr) = data(end,2);
            vis_Segment(cl+1,cl) = data(end,1);
        case char(9711) % Roundabout
            vis_Segment(2,3:end-1) = flipud(data(1:end-2,1));
            vis_Segment(3:end-1,end-1) = flipud(data(1:end-2,2));
            vis_Segment(end-1,2:end-2) = data(1:end-2,3);
            vis_Segment(2:end-2,2) = data(1:end-2,4);
            vis_Segment(1,cl) = data(end-1,1);
            vis_Segment(1,cr) = data(end,1);
            vis_Segment(cl,end) = data(end-1,2);
            vis_Segment(cr,end) = data(end,2);
            vis_Segment(end,cr) = data(end-1,3);
            vis_Segment(end,cl) = data(end,3);
            vis_Segment(cr,1) = data(end-1,4);
            vis_Segment(cl,1) = data(end,4);
        case char(9532)
            ls = segLength/2-3/2; % Length of the side roads
            tl = data(1:ls,1);
            tr = data(ls+1:end-1,1);
            rl = data(1:ls,2);
            rr = data(ls+1:end-1,2);
            bl = data(1:ls,3);
            br = data(ls+1:end-1,3);
            ll = data(1:ls,4);
            lr = data(ls+1:end-1,4);
            tc = data(end,1);
            rc = data(end,2);
            bc = data(end,3);
            lc = data(end,4);
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
%         case 'EX'
%             data = inOut_T{i,j}+0.2;
%             if i-1 == 0 % Top side
%             elseif ~strcmp(road_Network(i-1,j),'')
%                 vis_Segment(1,cr) = data(1,1);
%             end
%             if j+1 > simSize % Right side
%             elseif ~strcmp(road_Network(i,j+1),'')
%                 vis_Segment(cr,end) = data(1,2);
%             end
%             if i+1 > simSize % Bottom side
%             elseif ~strcmp(road_Network(i+1,j),'')
%                 vis_Segment(end,cl) = data(1,3);
%             end
%             if j-1 == 0 % Left side
%             elseif ~strcmp(road_Network(i,j-1),'')
%                 vis_Segment(cl,1) = data(1,4);
%             end
        otherwise % Case of empty cell ''
            
    end
                
    for i = 1:segLength
        for j = 1:segLength
            bID = i+(segLength-j)*segLength; % Button ID number
            tag = sprintf('edit%i',bID); % Generate tag for button
            [loci, locj] = quorem(sym(bID-1), sym(segLength)); % Calculate location in standard matrix notation based on button number
            loci = loci+1;
            locj = locj+1;
            vis = 'on';
            if isnan(vis_Segment(loci,locj))
                vis = 'off';
            end
            uicontrol('Parent',hObject,'Style','edit','String',vis_Segment(loci,locj),'Units','normalized',...
                'Position',[0.15+dist*(i-1) 0.25+dist*(j-1) boxSize boxSize],'Visible',vis,...
                'Callback',{@change_Value,hObject},'Tag', tag); % Create pushbutton
        end
    end
    
    uicontrol('Parent',hObject,'Style','pushbutton','String','Resume Simulation!',...
        'Units','normalized','Position',[0.3 0.05 0.4 0.15],'Visible','on',...
        'Callback',{@save_Changes,hObject},'Tag', 'Start'); % Create start simulation button

end
handles.Data = vis_Segment;
% Update handles structure
guidata(hObject, handles);


% UIWAIT makes CapandFlow_Edit wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = CapandFlow_Edit_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function change_Value(src,event,hObject)

handles = guidata(hObject);
segLength = handles.segmentLength;

num = src.Tag;
num = str2double(num(5:end));
[loci, locj] = quorem(sym(num-1), sym(segLength)); % Calculate location in standard matrix notation based on button number
loci = loci+1;
locj = locj+1;

data = handles.Data;
data(loci,locj) = str2double(src.String);

handles.Data = data;

guidata(hObject,handles);




function save_Changes(src,event,hObject)
handles = guidata(hObject);
segLength = handles.segmentLength;
flow = handles.Flow;
capacity = handles.Capacity;
Segi = handles.Segi;
Segj = handles.Segj;
CorF = handles.CorF;
road_Network = handles.Network;
data = handles.Data;

cl = segLength/2-0.5; % Constant, position of left lane
cr = segLength/2+1.5; % Constant, position of right lane

outData = zeros(size(flow{Segi,Segj}));

switch road_Network(Segi,Segj)
    case char(9472) % Straight L-R
        outData(:,1) = fliplr(data(cl,:));
        outData(:,2) = data(cr,:);
    case char(9474) % Straight T-B
        outData(:,1) = data(:,cl);
        outData(:,2) = flipud(data(:,cr));
    case char(9488) % Corner L-B
        outData(1:cl,2) = data(cr,1:cl);
        outData(cl+1:end,2) = data(cr+1:end,cl);
        outData(1:cl-2,1) = flipud(data(cr+2:end,cr));
        outData(cl-1,1) = flipud(data(cr,cr));
        outData(cl,1) = flipud(data(cl,cl+1));
        outData(cl+1,1) = fliplr(data(cl,cl));
        outData(cl+2:end,1) = fliplr(data(cl,1:cl-2));
    case char(9484) % Corner B-R
        outData(1:cl,2) = flipud(data(cr:end,cr));
        outData(cl+1:end,2) = data(cr,cr+1:end);
        outData(cl+2:end,1) = data(cr+2:end,cl);
        outData(cl+1,1) = data(cr,cl);
        outData(cl,1) = data(cl,cl+1);
        outData(cl-1,1) = data(cl,cr);
        outData(1:cl-2,1) = fliplr(data(cl,cr+2:end));
    case char(9492) % Corner R-T
        outData(1:cl,2) = fliplr(data(cl,cr:end));
        outData(cl+1:end,2) = flipud(data(1:cl-1,cr));
        outData(1:cl-2,1) = data(1:cl-2,cl);
        outData(cl-1,1) = data(cl,cl);
        outData(cl,1) = data(cr,cl+1);
        outData(cl+1,1) = data(cr,cr);
        outData(cl+2:end,1) = data(cr,cr+2:end);
    case char(9496) % Corner T-L
        outData(1:cl,2) = data(1:cl,cl);
        outData(cl+1:end,2) = fliplr(data(cl,1:cl-1));
        outData(1:cl-2,1) = data(cr,1:cl-2);
        outData(cl-1,1) = data(cr,cl);
        outData(cl,1) = data(cl+1,cr);
        outData(cl+1,1) = data(cl,cr);
        outData(cl+2:end,1) = flipud(data(1:cl-2,cr));
    case char(9508) % Intersection Left
        outData(1:segLength,1) = data(:,cl);
        outData(1:segLength,2) = flipud(data(:,cr));
        outData(segLength+1:end-1,2) = fliplr(data(cl,1:cl-1));
        outData(segLength+1:end-1,1) = data(cr,1:cl-1);
        outData(end,2) = data(cl,cl+1);
        outData(end,1) = data(cr,cl+1);
    case char(9516) % Intersection Bottom
        outData(1:segLength,2) = fliplr(data(cl,:));
        outData(1:segLength,1) = data(cr,:);
        outData(segLength+1:end-1,2) = data(cr+1:end,cl);
        outData(segLength+1:end-1,1) = flipud(data(cr+1:end,cr));
        outData(end,2) = data(cl+1,cl);
        outData(end,1) = data(cl+1,cr);
    case char(9500) % Intersection Right
        outData(1:segLength,2) = data(:,cl);
        outData(1:segLength,1) = flipud(data(:,cr));
        outData(segLength+1:end-1,2) = data(cr,cr+1:end);
        outData(segLength+1:end-1,1) = fliplr(data(cl,cr+1:end));
        outData(end,2) = data(cr,cl+1);
        outData(end,1) = data(cl,cl+1);
    case char(9524) % Intersection Top
        outData(1:segLength,2) = data(cr,:);
        outData(1:segLength,1) = fliplr(data(cl,:));
        outData(segLength+1:end-1,2) = flipud(data(1:cl-1,cr));
        outData(segLength+1:end-1,1) = data(1:cl-1,cl);
        outData(end,2) = data(cl+1,cr);
        outData(end,1) = data(cl+1,cl);
    case char(9711) % Roundabout
        outData(1:end-2,1) = fliplr(data(2,3:end-1));
        outData(1:end-2,2) = flipud(data(3:end-1,end-1));
        outData(1:end-2,3) = data(end-1,2:end-2);
        outData(1:end-2,4) = data(2:end-2,2);
        outData(end-1,1) = data(1,cl);
        outData(end,1) = data(1,cr);
        outData(end-1,2) = data(cl,end);
        outData(end,2) = data(cr,end);
        outData(end-1,3) = data(end,cr);
        outData(end,3) = data(end,cl);
        outData(end-1,4) = data(cr,1);
        outData(end,4) = data(cl,1);
    case char(9532) % Traffic Light % NOT DONE!
        ls = segLength/2-3/2; % Length of the side roads
        tl = data(1:ls,cl); %Done
        tr = data(1:ls+1,cr); %Done
        rl = data(cl,end-(ls-1):end); %Done
        rr = data(cr,end-ls:end); %Done
        bl = data(end-(ls-1):end,cr); %Done
        br = data(end-ls:end,cl); %Done
        ll = data(cr,1:ls); %Done
        lr = data(cl,1:ls+1); %Done
        tc = data(cl,cl+1); %Done
        rc = data(cl+1,cr); %Done
        bc = data(cr,cl+1); %Done
        lc = data(cl+1,cl); %Done
        outData(1:ls,1) = tl; %Done
        outData(1:ls,2) = fliplr(rl); %Done
        outData(1:ls,3) = flipud(bl); %Done
        outData(1:ls,4) = ll; %Done
        outData(ls+1:end-1,1) = flipud(tr); %Done
        outData(ls+1:end-1,2) = rr; %Done
        outData(ls+1:end-1,3) = br; %Done
        outData(ls+1:end-1,4) = fliplr(lr); %Done
        outData(end,1) = tc; %Done
        outData(end,2) = rc; %Done
        outData(end,3) = bc; %Done
        outData(end,4) = lc; %Done
%     case 'EX'
%         data = inOut_T{i,j}+0.2;
%         if i-1 == 0 % Top side
%         elseif ~strcmp(road_Network(i-1,j),'')
%             vis_Segment(1,cr) = data(1,1);
%         end
%         if j+1 > simSize % Right side
%         elseif ~strcmp(road_Network(i,j+1),'')
%             vis_Segment(cr,end) = data(1,2);
%         end
%         if i+1 > simSize % Bottom side
%         elseif ~strcmp(road_Network(i+1,j),'')
%             vis_Segment(end,cl) = data(1,3);
%         end
%         if j-1 == 0 % Left side
%         elseif ~strcmp(road_Network(i,j-1),'')
%             vis_Segment(cl,1) = data(1,4);
%         end
    otherwise % Case of empty cell ''

end

if strcmp(CorF,'Capacity')
    capacity(Segi,Segj) = {outData};
    handles.Capacity = capacity;
else
    flow(Segi,Segj) = {outData};
    handles.Flow = flow;
end

CPWindow = handles.CPWindow;
CPHandles = guidata(CPWindow);
CPHandles.Network = handles.Network;
CPHandles.Connections = handles.Connections;
CPHandles.Segments = handles.Segments;
CPHandles.Capacity = handles.Capacity;
CPHandles.Flow = handles.Flow;
guidata(CPWindow,CPHandles);
close;
