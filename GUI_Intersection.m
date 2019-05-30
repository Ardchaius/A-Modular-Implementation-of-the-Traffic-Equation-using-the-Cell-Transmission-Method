% Improved Intersection Segment function
% IMPROVED AND CORRECT OVER GUI_Intersection()!!!

function [mr1,ml1,sr1,sl1,t11,t21] = GUI_Intersection(mr,ml,sr,sl,t1,t2,lm)
% mr- Main road right lane including input and output cells
% ml- Main road left lane including input and output cells
% sr- Side road right lane including only output cell
% sl- Side road left lane including only input cell
% t1- Turn off from mr to sr no input or output cells
% t2- Turn off from sl to mr no input or output cells
% lm- length of main road

%% Constants

ls = lm/2-3/2; % Length of side road based on length of main road
v = 1; % Forward Velocity
w = 0.75; % Backward wave velocity
nN = 3; % Max vehicles in each cell
qQ = 3; % Maximum flow through each cell
ill = lm/2 + 0.5; % Intersection location of left side road (adjusted for input and output cells) with respect to main right
ilr = lm/2 + 2.5; % Intersection location of right side road (adjusted for input and output cells)
sl2mr = 0.5; % Percentage of cars going from sl to mr
mr2sr = 0.5; % Percentage of cars going from mr to sr
ml2sr = 0.5; % Percentage of cars going from ml to sr


%% Initilization

Nm = nN * ones(lm+2,1); % Max vehicles for each cell main road
Ns = nN * ones(ls+2,1); % Max vehicles for each cell side road

Qm = qQ * ones(lm+2,1); % Vectorized Q for use in min functions for main road
Qs = qQ * ones(ls+2,1); % Vectorized Q for use in min functions for side road

sr = [ml(ill)*ml2sr;sr]; % Left and not right (ill not ilr) due  to the direction of travel
ml(ill) = ml(ill) - sr(1);

t1 = [mr2sr*mr(ill+1);t1;sr(2)];
mr(ill+1) = mr(ill+1) - t1(1);

sl = [sl;ml(ilr)];
t2 = [sl2mr*sl(end-2);t2;mr(ill)];
sl(end-2) = sl(end-2) - t2(1);

%% Calculations

Smr = min([Qm,mr],[],2); % Calculate the send capacity for each cell right lane of main road
Rmr = min([Qm,(w/v)*(Nm-mr)],[],2); % Calculate the recieve capacity for each cell right lane of main road

Sml = min([Qm,ml],[],2); % Calculate the send capacity for each cell left lane of main road
Rml = min([Qm,(w/v)*(Nm-ml)],[],2); % Calculate the recieve capacity for each cell left lane of main road

Ssr = min([Qs,sr],[],2); % Calculate the send capacity for each cell right lane of side road
Rsr = min([Qs,(w/v)*(Ns-sr)],[],2); % Calculate the recieve capacity for each cell right lane of side road

Ssl = min([Qs,sl],[],2); % Calculate the send capacity for each cell left lane of side road
Rsl = min([Qs,(w/v)*(Ns-sl)],[],2); % Calculate the recieve capacity for each cell left lane of side road

ymr = min([Smr(1:end-1),Rmr(2:end)],[],2); % Calculate flow into each cell from cell 2 to end main right
yml = min([Sml(1:end-1),Rml(2:end)],[],2); % Calculate flow into each cell from cell 2 to end main left
ysr = min([Ssr(1:end-1),Rsr(2:end)],[],2); % Calculate flow into each cell from cell 2 to end side right
ysl = min([Ssl(1:end-1),Rsl(2:end)],[],2); % Calculate flow into each cell from cell 2 to end side left

St1 = min([[qQ;1;qQ],t1],[],2);
Rt1 = min([[qQ;1;qQ],(w/v)*([nN;1;nN]-t1)],[],2);

yt1 = min([St1(1:end-1),Rt1(2:end)],[],2); % Calculate flow for t1 section
yt1(end) = min([yt1(end),nN - ml(ill) - yml(ill-1)]);

St2 = min([[qQ;1;qQ],t2],[],2);
Rt2 = min([[qQ;1;qQ],(w/v)*([nN;1;nN]-t2)],[],2);

yt2 = min([St2(1:end-1),Rt2(2:end)],[],2); % Calculate flow for t2 section
yt2(end) = min([yt2(end),nN - ml(ilr) - yml(ilr-1),nN - mr(ill) - ymr(ill-1)]);

mr1 = mr(2:end-1) + ymr(1:end-1) - ymr(2:end);
ml1 = ml(2:end-1) + yml(1:end-1) - yml(2:end);
sl1 = sl(2:end-1) + ysl(1:end-1) - ysl(2:end);

sr1 = sr;
sr1(1:end-1) = sr1(1:end-1) - ysr;
sr1(2:end) = sr1(2:end) + ysr;

t21 = t2;
t21(1:2) = t21(1:2) - yt2;
t21(2:3) = t21(2:3) + yt2;

t11 = t1;
t11(1:2) = t11(1:2) - yt1;
t11(2:3) = t11(2:3) + yt1;

sl1(end-1) = sl1(end-1) + t21(1);
mr1(ill-1) = mr1(ill-1) + yt2(end);
mr1(ill+1) = mr1(ill+1) + t11(1);
sr1(2) = sr1(2) + yt1(end);
ml1(ill-1) = ml1(ill-1) + sr1(1);
ml1(ill+1) = ml1(ill+1) + ysl(end);

sr1 = sr1(2:end-1);
t11 = t11(2);
t21 = t21(2);

end