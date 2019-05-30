function [tl_1,tr_1,rl_1,rr_1,bl_1,br_1,ll_1,lr_1,tc_1,rc_1,bc_1,lc_1,timer] = GUI_Traffic_Light(tl,tr,rl,rr,bl,br,ll,lr,tc,rc,bc,lc,length,timer)

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
lm = length; % Length of main road
ls = length/2-3/2; % Length of the side roads
v = 1; % Forward Velocity
w = 0.75; % Backward wave velocity
nN = 3; % Max vehicles in each cell
qQ = 3; % Maximum flow through each cell
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

Nm = nN * ones(lm+2,1); % Max vehicles for each cell main road
Ns = nN * ones(ls+2,1); % Max vehicles for each cell side road
Nc = nN * ones(3,1); % Max vehicles for each cell in crossing

Qm = qQ * ones(lm+2,1); % Vectorized Q for use in min functions for main road
Qsr = qQ * ones(ls+2,1); % Vectorized Q for use in min functions for side road right side
Qsl = qQ * ones(ls+2,1); % Vectorized Q for use in min functions for side road left side
Qsl(end) = 0;
Qc = qQ * ones(3,1); % Vectorized Q for use in min functions for each cell in crossing

if (0 < timer) && (timer < light_Switch + light_Pause + 1) % Up down lanes have green light
    mr = [bl;rr(1);rc;tr];
    ml = [tl;lr(1);lc;br];
    slr = lr(2:end); % Output to the left, seen from the perspective of the up lane
    srr = rr(2:end); % Output to the right, seen from the perspective of the up lane
    sll = [ll;br(1)]; % Input from the left, seen from the perspective of the up lane
    srl = [rl;tr(1)]; % Input from the right, seen from the perspective of the up lane
    topc = tc;
    botc = bc;
    mrsltr = u2l; % Turning rate from main right to side left road
    mrsrtr = u2r; % Turning rate from main right to side right road
    mlsltr = d2l; % Turning rate from main left to side left road
    mlsrtr = d2r; % Turning rate from main left to side right road
    if (light_Switch < timer)
        Qm(ls+2) = 0;
    end
    [mr_1,ml_1,slr_1,sll_1,srr_1,srl_1,topc_1,botc_1] = GUI_Traffic_Light_Sub(mr,ml,slr,sll,srr,srl,topc,botc,mrsltr,mrsrtr,mlsltr,mlsrtr,ls,v,w,nN,Nm,Ns,Nc,Qm,Qsl,Qsr,Qc);
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
    mrsltr = r2l; % Turning rate from main right to side left road
    mrsrtr = r2r; % Turning rate from main right to side right road
    mlsltr = l2l; % Turning rate from main left to side left road
    mlsrtr = l2r; % Turning rate from main left to side right road
    if (2 * light_Switch + light_Pause < timer)
        Qm(ls+2) = 0;
    end
    [mr_1,ml_1,slr_1,sll_1,srr_1,srl_1,topc_1,botc_1] = GUI_Traffic_Light_Sub(mr,ml,slr,sll,srr,srl,topc,botc,mrsltr,mrsrtr,mlsltr,mlsrtr,ls,v,w,nN,Nm,Ns,Nc,Qm,Qsl,Qsr,Qc);
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

function [mr_1,ml_1,slr_1,sll_1,srr_1,srl_1,topc_1,botc_1] = GUI_Traffic_Light_Sub(mr,ml,slr,sll,srr,srl,topc,botc,mrsltr,mrsrtr,mlsltr,mlsrtr,ls,v,w,nN,Nm,Ns,Nc,Qm,Qsl,Qsr,Qc)

    %% Preparations
    slr = [mlsltr*ml(ls+2);slr];
    ml(ls+2) = ml(ls+2) - slr(1);
    srr = [mrsrtr*mr(ls+2);srr];
    mr(ls+2) = mr(ls+2) - srr(1);
    
    topc = [mrsltr*mr(ls+4);topc;slr(2)];
    mr(ls+4) = mr(ls+4) - topc(1);
    botc = [mlsrtr*ml(ls+4);botc;srr(2)];
    ml(ls+4) = ml(ls+4) - botc(1);
    
    
    %% Calculations
    Sml = min([Qm,ml],[],2); % Calculate the send capacity for each cell left lane of main road
    Rml = min([Qm,(w/v)*(Nm-ml)],[],2); % Calculate the recieve capacity for each cell left lane of main road
    Smr = min([Qm,mr],[],2); % Calculate the send capacity for each cell right lane of main road
    Rmr = min([Qm,(w/v)*(Nm-mr)],[],2); % Calculate the recieve capacity for each cell right lane of main road
    
    Ssll = min([Qsl,sll],[],2);
    Rsll = min([Qsl,(w/v)*(Ns-sll)],[],2);
    Ssrl = min([Qsl,srl],[],2);
    Rsrl = min([Qsl,(w/v)*(Ns-srl)],[],2);
    
    Sslr = min([Qsr,slr],[],2);
    Rslr = min([Qsr,(w/v)*(Ns-slr)],[],2);
    Ssrr = min([Qsr,srr],[],2);
    Rsrr = min([Qsr,(w/v)*(Ns-srr)],[],2);
    
    yml = min([Sml(1:end-1),Rml(2:end)],[],2);
    ymr = min([Smr(1:end-1),Rmr(2:end)],[],2);
    ysll = min([Ssll(1:end-1),Rsll(2:end)],[],2);
    ysrl = min([Ssrl(1:end-1),Rsrl(2:end)],[],2);
    yslr = min([Sslr(1:end-1),Rslr(2:end)],[],2);
    ysrr = min([Ssrr(1:end-1),Rsrr(2:end)],[],2);
    
    Stopc = min([Qc,topc],[],2);
    Rtopc = min([Qc,(w/v)*(Nc-topc)],[],2);
    Sbotc = min([Qc,botc],[],2);
    Rbotc = min([Qc,(w/v)*(Nc-botc)],[],2);
    
    ytopc = min([Stopc(1:end-1),Rtopc(2:end)],[],2);
    ytopc(end) = min([ytopc(end),nN - slr(2) - yslr(1)],[],2);
    
    ybotc = min([Sbotc(1:end-1),Rbotc(2:end)],[],2);
    ybotc(end) = min([ybotc(end),nN - srr(2) - ysrr(1)],[],2);
    
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