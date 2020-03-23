    % Improved Straight Segment function

function [segmentl_t1, segmentr_t1] = GUI_Straight_CF(segmentl_t, segmentr_t, mlt_C, mrt_C, mlt_F, mrt_F)

%% Calculations

    Sl = min([mlt_F,segmentl_t],[],2); % Calculate the send capacity for each cell
    Rl = max(min([mlt_F,(mlt_C-segmentl_t)],[],2),0); % Calculate the recieve capacity for each cell

    yl = min([Sl(1:end-1),Rl(2:end)],[],2); % Calculate flow into each cell from cell 2 to end

    segmentl_t1 = segmentl_t(2:end-1) + yl(1:end-1) - yl(2:end); % Calculate next time step
    
    Sr = min([mrt_F,segmentr_t],[],2); % Calculate the send capacity for each cell
    Rr = max(min([mrt_F,(mrt_C-segmentr_t)],[],2),0); % Calculate the recieve capacity for each cell

    yr = min([Sr(1:end-1),Rr(2:end)],[],2); % Calculate flow into each cell from cell 2 to end

    segmentr_t1 = segmentr_t(2:end-1) + yr(1:end-1) - yr(2:end); % Calculate next time step

end