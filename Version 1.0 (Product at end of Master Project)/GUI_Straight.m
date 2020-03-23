% Improved Straight Segment function

function [segmentl_t1, segmentr_t1] = GUI_Straight(segmentl_t, segmentr_t, l)

%% Constants

v = 1; % Forward Velocity
w = 0.75; % Backward wave velocity
nN = 3; % Max vehicles in each cell
qQ = 3; % Maximum flow through each cell


%% Initializations

N = nN * ones(l+2,1); % Max vehicles for each cell

Q = qQ * ones(l+2,1); % Vectorized Q for use in min functions



%% Calculations

    Sl = min([Q,segmentl_t],[],2); % Calculate the send capacity for each cell
    Rl = min([Q,(w/v)*(N-segmentl_t)],[],2); % Calculate the recieve capacity for each cell

    yl = min([Sl(1:end-1),Rl(2:end)],[],2); % Calculate flow into each cell from cell 2 to end

    segmentl_t1 = segmentl_t(2:end-1) + yl(1:end-1) - yl(2:end); % Calculate next time step
    
    Sr = min([Q,segmentr_t],[],2); % Calculate the send capacity for each cell
    Rr = min([Q,(w/v)*(N-segmentr_t)],[],2); % Calculate the recieve capacity for each cell

    yr = min([Sr(1:end-1),Rr(2:end)],[],2); % Calculate flow into each cell from cell 2 to end

    segmentr_t1 = segmentr_t(2:end-1) + yr(1:end-1) - yr(2:end); % Calculate next time step

end