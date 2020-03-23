function [tl_1,tr_1,rl_1,rr_1,bl_1,br_1,ll_1,lr_1,tc_1,rc_1,bc_1,lc_1,timer] = GUI_Traffic_Light_CF...
    (tl,tr,rl,rr,bl,br,ll,lr,tc,rc,bc,lc,...
    tl_C,tr_C,rl_C,rr_C,bl_C,br_C,ll_C,lr_C,tc_C,rc_C,bc_C,lc_C,...
    tl_F,tr_F,rl_F,rr_F,bl_F,br_F,ll_F,lr_F,tc_F,rc_F,bc_F,lc_F,length,timer)

% tl- Top side left lane, includes input from top
% tr- Top side right lane plus cell in the intersection itself, includes
% output to top
% rl- Right side left lane, includes input from right
% rr- Right side right lane plus cell in the intersection itself, includes
% output to right
% bl- Bottom side left lane, includes input from bottom
% br- Bottom side right lane plus cell in the intersection itself, includes
% output to bottom
% ll- Left side left lane, includes input from left
% lr- Left side right lane plus cell in the intersection itself, includes
% output to left
% length- Length on the segment
% timer- Value to keep track of which traffic lights are red and for how
% long

%% Constants
ls = length/2-3/2; % Length of the side roads
v = 1; % Forward Velocity
w = 0.75; % Backward wave velocity
u2l = 1/3; % Turning rate from up lane to left lane
u2r = 1/4; % Turning rate from up lane to right lane 
d2l = 1/3; % Turning rate from down lane to right lane
d2r = 1/4; % Turning rate from down lane to left lane
r2r = 1/4; % Turning rate from right lane to down lane
r2l = 1/3; % Turning rate from right lane to up lane
l2l = 1/4; % Turning rate from left lane to up lane
l2r = 1/3; % Turning rate from left lane to down lane

light_Switch = 10; % Number of iterations before light changes
light_Pause = 3; % Number of iterations were both lights are red, clearing the 


%% Initializations
% Here we find out what light is currently green and generalize the code
% accordingly. This neatens the code, although it might appear confusing.
% The situation is the same, no matter which direction currently has the
% green light.

if (0 < timer) && (timer < light_Switch + light_Pause + 1) % Up down lanes have green light
    mr = [bl;rr(1);rc;tr];
    ml = [tl;lr(1);lc;br];
    slr = lr(2:end); % Output to the left, seen from the perspective of the up lane
    srr = rr(2:end); % Output to the right, seen from the perspective of the up lane
    sll = [ll;br(1)]; % Input from the left, seen from the perspective of the up lane
    srl = [rl;tr(1)]; % Input from the right, seen from the perspective of the up lane
    topc = tc;
    botc = bc;
    mr_C = [bl_C;rr_C(1);rc_C;tr_C];
    ml_C = [tl_C;lr_C(1);lc_C;br_C];
    slr_C = lr_C(2:end);
    srr_C = rr_C(2:end);
    sll_C = [ll_C;br_C(1)];
    srl_C = [rl_C;tr_C(1)];
    topc_C = tc_C;
    botc_C = bc_C;
    mr_F = [bl_F;rr_F(1);rc_F;tr_F];
    ml_F = [tl_F;lr_F(1);lc_F;br_F];
    slr_F = lr_F(2:end);
    srr_F = rr_F(2:end); 
    sll_F = [ll_F;0]; 
    srl_F = [rl_F;0];
    topc_F = tc_F;
    botc_F = bc_F;
    
    mrsltr = u2l; % Turning rate from main right to side left road
    mrsrtr = u2r; % Turning rate from main right to side right road
    mlsltr = d2l; % Turning rate from main left to side left road
    mlsrtr = d2r; % Turning rate from main left to side right road
    if (light_Switch < timer)
        mr_F(ls+1) = 0;
        ml_F(ls+1) = 0;
    end
    % Calculate next time step
    [mr_1,ml_1,slr_1,sll_1,srr_1,srl_1,topc_1,botc_1] = GUI_Traffic_Light_Sub_CF...
        (mr,ml,slr,sll,srr,srl,topc,botc,...
        mr_C,ml_C,slr_C,sll_C,srr_C,srl_C,topc_C,botc_C,...
        mr_F,ml_F,slr_F,sll_F,srr_F,srl_F,topc_F,botc_F,...
        mrsltr,mrsrtr,mlsltr,mlsrtr,ls,v,w);
    % Reformat data for output
    tl_1 = ml_1(1:ls);
    tr_1 = mr_1(end-ls:end);
    rl_1 = srl_1;
    rr_1 = [mr_1(ls+1);srr_1];
    bl_1 = mr_1(1:ls);
    br_1 = ml_1(end-ls:end);
    ll_1 = sll_1;
    lr_1 = [ml_1(ls+1);slr_1];
    tc_1 = topc_1;
    rc_1 = mr_1(ls+2);
    bc_1 = botc_1;
    lc_1 = ml_1(ls+2);
    timer = timer + 1;
elseif (light_Switch + light_Pause < timer) && (timer < 2 * (light_Switch + light_Pause) + 1) % Right left lanes have green light    
    mr = [ll;br(1);bc;rr];
    ml = [rl;tr(1);tc;lr];
    slr = tr(2:end); % Output to the left, seen from the perspective of the right lane
    srr = br(2:end); % Output to the right, seen from the perspective of the up lane
    sll = [tl;lr(1)]; % Input from the left, seen from the perspective of the right lane
    srl = [bl;rr(1)]; % Input from the right, seen from the perspective of the right lane
    topc = rc;
    botc = lc;
    mr_C = [ll_C;br_C(1);bc_C;rr_C];
    ml_C = [rl_C;tr_C(1);tc_C;lr_C];
    slr_C = tr_C(2:end); 
    srr_C = br_C(2:end); 
    sll_C = [tl_C;lr_C(1)]; 
    srl_C = [bl_C;rr_C(1)]; 
    topc_C = rc_C;
    botc_C = lc_C;
    mr_F = [ll_F;br_F(1);bc_F;rr_F];
    ml_F = [rl_F;tr_F(1);tc_F;lr_F];
    slr_F = tr_F(2:end); 
    srr_F = br_F(2:end); 
    sll_F = [tl_F;0];
    srl_F = [bl_F;0];
    topc_F = rc_F;
    botc_F = lc_F;
    mrsltr = r2l; % Turning rate from main right to side left road
    mrsrtr = r2r; % Turning rate from main right to side right road
    mlsltr = l2l; % Turning rate from main left to side left road
    mlsrtr = l2r; % Turning rate from main left to side right road
    if (2 * light_Switch + light_Pause < timer)
        mr_F(ls+1) = 0;
        ml_F(ls+1) = 0;
    end
    % Calcualte next time step
    [mr_1,ml_1,slr_1,sll_1,srr_1,srl_1,topc_1,botc_1] = GUI_Traffic_Light_Sub_CF...
        (mr,ml,slr,sll,srr,srl,topc,botc,...
        mr_C,ml_C,slr_C,sll_C,srr_C,srl_C,topc_C,botc_C,...
        mr_F,ml_F,slr_F,sll_F,srr_F,srl_F,topc_F,botc_F,...
        mrsltr,mrsrtr,mlsltr,mlsrtr,ls,v,w);
    % Reformat data for output
    tl_1 = sll_1;
    tr_1 = [ml_1(ls+1);slr_1];
    rl_1 = ml_1(1:ls);
    rr_1 = mr_1(end-ls:end);
    bl_1 = srl_1;
    br_1 = [mr_1(ls+1);srr_1];
    ll_1 = mr_1(1:ls);
    lr_1 = ml_1(end-ls:end);
    tc_1 = ml_1(ls+2);
    rc_1 = topc_1;
    bc_1 = mr_1(ls+2);
    lc_1 = botc_1;
    timer = timer + 1;
    if timer > 2 * (light_Switch + light_Pause)
        timer = 1;
    end
end

end

function [mr_1,ml_1,slr_1,sll_1,srr_1,srl_1,topc_1,botc_1] = GUI_Traffic_Light_Sub_CF...
    (mr,ml,slr,sll,srr,srl,topc,botc,...
        mr_C,ml_C,slr_C,sll_C,srr_C,srl_C,topc_C,botc_C,...
        mr_F,ml_F,slr_F,sll_F,srr_F,srl_F,topc_F,botc_F,...
        mrsltr,mrsrtr,mlsltr,mlsrtr,ls,v,w)

    %% Preparations
    
    slr = [mlsltr*ml(ls+2);slr];
    slr_C = [ml_C(ls+2);slr_C];
    slr_F = [ml_F(ls+2);slr_F];
    ml(ls+2) = ml(ls+2) - slr(1);
    srr = [mrsrtr*mr(ls+2);srr];
    srr_C = [mr_C(ls+2);srr_C];
    srr_F = [mr_F(ls+2);srr_F];
    mr(ls+2) = mr(ls+2) - srr(1);
    
    topc = [mrsltr*mr(ls+4);topc;slr(2)];
    topc_C = [mr_C(ls+4);topc_C;slr_C(2)];
    topc_F = [mr_F(ls+4);topc_F;slr_F(2)];
    mr(ls+4) = mr(ls+4) - topc(1);
    botc = [mlsrtr*ml(ls+4);botc;srr(2)];
    botc_C = [ml_C(ls+4);botc_C;srr_C(2)];
    botc_F = [ml_F(ls+4);botc_F;srr_F(2)];
    ml(ls+4) = ml(ls+4) - botc(1);
    
    
    %% Calculations
    Sml = min([ml_F,ml],[],2); % Calculate the send capacity for each cell left lane of main road
    Rml = max(min([ml_F,(w/v)*(ml_C-ml)],[],2),0); % Calculate the recieve capacity for each cell left lane of main road
    Smr = min([mr_F,mr],[],2); % Calculate the send capacity for each cell right lane of main road
    Rmr = max(min([mr_F,(w/v)*(mr_C-mr)],[],2),0); % Calculate the recieve capacity for each cell right lane of main road
    
    Ssll = min([sll_F,sll],[],2); % Calculate send capacity for the left lane of the left side road
    Rsll = max(min([sll_F,(w/v)*(sll_C-sll)],[],2),0); % Calculate recieve capacity for the left lane of the left side road
    Ssrl = min([srl_F,srl],[],2); % Calculate send capacity for the left lane of the right side road
    Rsrl = max(min([srl_F,(w/v)*(srl_C-srl)],[],2),0); % Calculate recieve capacity for the left lane of the right side road
    
    Sslr = min([slr_F,slr],[],2); % Calculate send capacity for the right lane of the left side road
    Rslr = max(min([slr_F,(w/v)*(slr_C-slr)],[],2),0); % Calculate recieve capacity for the right lane of the left side road
    Ssrr = min([srr_F,srr],[],2); % Calcualte send capacity for the right lane of the right side road
    Rsrr = max(min([srr_F,(w/v)*(srr_C-srr)],[],2),0); % Calcualte recieve capacity for the right lane of the right side road
    
    yml = min([Sml(1:end-1),Rml(2:end)],[],2); % Calculate the flow of vehicles for all lanes
    ymr = min([Smr(1:end-1),Rmr(2:end)],[],2);
    ysll = min([Ssll(1:end-1),Rsll(2:end)],[],2);
    ysrl = min([Ssrl(1:end-1),Rsrl(2:end)],[],2);
    yslr = min([Sslr(1:end-1),Rslr(2:end)],[],2);
    ysrr = min([Ssrr(1:end-1),Rsrr(2:end)],[],2);
    
    Stopc = min([topc_F,topc],[],2); % Calculate send capacity for cell between the main lanes top
    Rtopc = max(min([topc_F,(w/v)*(topc_C-topc)],[],2),0); % Calculate recieve capacity for cell between the main lanes top
    Sbotc = min([botc_F,botc],[],2); % Calculate send capacity for cell between the main lanes bottom
    Rbotc = max(min([botc_F,(w/v)*(botc_C-botc)],[],2),0); % Calculate recieve capacity for cell between the main lanes bottom
    
    ytopc = min([Stopc(1:end-1),Rtopc(2:end)],[],2); % Calculate flow of vehicles for cell between main lanes top
    ytopc(end) = min([ytopc(end),slr_C(2) - slr(2) - yslr(1)],[],2);
    
    ybotc = min([Sbotc(1:end-1),Rbotc(2:end)],[],2); % Calculate flow of vehicles for cell between main lanes bottom
    ybotc(end) = min([ybotc(end),srr_C(2) - srr(2) - ysrr(1)],[],2);
    
    
    % Calculate new densities based on the flow of vehicles for all lanes
    ml_1 = ml(2:end-1) + yml(1:end-1) - yml(2:end);
    mr_1 = mr(2:end-1) + ymr(1:end-1) - ymr(2:end);
    
    sll_1 = sll(2:end-1) + ysll(1:end-1) - ysll(2:end);
    srl_1 = srl(2:end-1) + ysrl(1:end-1) - ysrl(2:end);
    
    slr_1 = slr;
    slr_1(2:end) = slr_1(2:end) + yslr;
    slr_1(1:end-1) = slr_1(1:end-1) - yslr;
    
    srr_1 = srr;
    srr_1(2:end) = srr_1(2:end) + ysrr;
    srr_1(1:end-1) = srr_1(1:end-1) - ysrr;
    
    topc_1 = topc;
    topc_1(2:end) = topc_1(2:end) + ytopc;
    topc_1(1:end-1) = topc_1(1:end-1) - ytopc;
    
    botc_1 = botc;
    botc_1(2:end) = botc_1(2:end) + ybotc;
    botc_1(1:end-1) = botc_1(1:end-1) - ybotc;
    
    ml_1(ls+1) = ml_1(ls+1) + slr_1(1);
    mr_1(ls+1) = mr_1(ls+1) + srr_1(1);
    ml_1(ls+3) = ml_1(ls+3) + botc_1(1);
    mr_1(ls+3) = mr_1(ls+3) + topc_1(1);
    
    slr_1 = slr_1(2:end-1);
    slr_1(1) = slr_1(1) + ytopc(end);
    srr_1 = srr_1(2:end-1);
    srr_1(1) = srr_1(1) + ybotc(end);
    topc_1 = topc_1(2:end-1);
    botc_1 = botc_1(2:end-1);

end