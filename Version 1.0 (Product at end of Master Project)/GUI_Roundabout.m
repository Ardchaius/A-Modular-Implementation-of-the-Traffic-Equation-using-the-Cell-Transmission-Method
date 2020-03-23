% Improved Roundabout Segment function

function [mt_1,mr_1,mb_1,ml_1,stl_1,str_1,srl_1,srr_1,sbl_1,sbr_1,sll_1,slr_1]...
    = GUI_Roundabout(mt,mr,mb,ml,stl,str,srl,srr,sbl,sbr,sll,slr,length)

%% Constants
    v = 1;
    w = 0.75; % Backward wave velocity
    nN = 3; % Max vehicles in each cell
    qQ = 3; % Maximum flow through each cell
    m_length = length-3; % Length of main road pieces
    s_length = 1; % Length of side road pieces
    exl = length/2 - 0.5; % Location of exit lane in direction of travel adjusted for input and output cells
    enl = length/2 + 1.5; % Location of entrance lane in direction of travel adjusted for input and output cells
    ert = 0.25; % Exit rate of top piece
    err = 0.25; % Exit rate of right piece
    erb = 0.25; % Exit rate of bottom piece
    erl = 0.25; % Exit rate of left piece
    
%% Initializations
    Nm = nN * ones(m_length+2,1); % Max vehicles for each cell main road
    Ns = nN * ones(s_length+2,1); % Max vehicles for each cell side road
    
    Qm = qQ * ones(m_length+2,1); % Vectorized Q for use in min functions for main road
    Qs = qQ * ones(s_length+2,1); % Vectorized Q for use in min functions for side road
    
    
    % Format data to correct format for use in calculations
    stl = [stl;mt(enl)];
    str = [mt(exl)*ert;str];
    mt(exl) = mt(exl)*(1-ert);
    
    srl = [srl;mr(enl)];
    srr = [mr(exl)*err;srr];
    mr(exl) = mr(exl)*(1-err);
    
    sbl = [sbl;mb(enl)];
    sbr = [mb(exl)*erb;sbr];
    mb(exl) = mb(exl)*(1-erb);
    
    sll = [sll;ml(enl)];
    slr = [ml(exl)*erl;slr];
    ml(exl) = ml(exl)*(1-erl);
    
%% Calculations
    [mt_1,stl_1,str_1] = GUI_Roundabout_Sub(mt,stl,str,enl,exl,w,v,Qm,Nm,Qs,Ns);
    [mr_1,srl_1,srr_1] = GUI_Roundabout_Sub(mr,srl,srr,enl,exl,w,v,Qm,Nm,Qs,Ns);
    [mb_1,sbl_1,sbr_1] = GUI_Roundabout_Sub(mb,sbl,sbr,enl,exl,w,v,Qm,Nm,Qs,Ns);
    [ml_1,sll_1,slr_1] = GUI_Roundabout_Sub(ml,sll,slr,enl,exl,w,v,Qm,Nm,Qs,Ns);

end

function [m_t1,sl_t1,sr_t1] = GUI_Roundabout_Sub(m,sl,sr,enl,exl,w,v,Qm,Nm,Qs,Ns)

    Sm = min([Qm,m],[],2); % Send and recieve main road top piece
    Rm = min([Qm,(w/v)*(Nm-m)],[],2);
    
    ym = min([Sm(1:end-1),Rm(2:end)],[],2);
    
    sl(end) = sl(end) + ym(enl-1);
    Ssl = min([Qs,sl],[],2); % Send and recieve for left and right lanes of the incoming and outgoing side roads
    Rsl = min([Qs,(w/v)*(Ns-sl)],[],2);
    Ssr = min([Qs,sr],[],2);
    Rsr = min([Qs,(w/v)*(Ns-sr)],[],2);
    
    ysl = min([Ssl(1:end-1),Rsl(2:end)],[],2); % Calculate veicle flow
    ysr = min([Ssr(1:end-1),Rsr(2:end)],[],2);
    
    % Calculate new time step
    m_t1 = m(2:end-1) + ym(1:end-1) - ym(2:end);
    m_t1(exl-1) = m_t1(exl-1) + sr(1) - ysr(1);
    m_t1(enl-1) = m_t1(enl-1) + ysl(end);
    
    sl_t1 = sl(2:end-1) + ysl(1:end-1) - ysl(2:end);
    sr_t1 = sr(2:end-1) + ysr(1:end-1) - ysr(2:end);

end