% Improved Roundabout Segment function

function [mt_1,mr_1,mb_1,ml_1,stl_1,str_1,srl_1,srr_1,sbl_1,sbr_1,sll_1,slr_1]...
    = GUI_Roundabout_CF(mt,mr,mb,ml,stl,str,srl,srr,sbl,sbr,sll,slr,...
                    mt_C,mr_C,mb_C,ml_C,stl_C,str_C,srl_C,srr_C,sbl_C,sbr_C,sll_C,slr_C,...
                    mt_F,mr_F,mb_F,ml_F,stl_F,str_F,srl_F,srr_F,sbl_F,sbr_F,sll_F,slr_F,length)

%% Constants
    v = 1;
    w = 0.75; % Backward wave velocity
    exl = length/2 - 0.5; % Location of exit lane in direction of travel adjusted for input and output cells
    enl = length/2 + 1.5; % Location of entrance lane in direction of travel adjusted for input and output cells
    ert = 0.25; % Exit rate of top piece
    err = 0.25; % Exit rate of right piece
    erb = 0.25; % Exit rate of bottom piece
    erl = 0.25; % Exit rate of left piece
    
%% Initializations
    
    % Format data to correct format for use in calculations
    stl = [stl;mt(enl)];
    str = [mt(exl)*ert;str];
    mt(exl) = mt(exl)*(1-ert);
    stl_C = [stl_C;mt_C(enl)];
    str_C = [mt_C(exl);str_C];
    stl_F = [stl_F;mt_F(enl)];
    str_F = [mt_F(exl);str_F];
    
    srl = [srl;mr(enl)];
    srr = [mr(exl)*err;srr];
    mr(exl) = mr(exl)*(1-err);
    srl_C = [srl_C;mr_C(enl)];
    srr_C = [mr_C(exl);srr_C];
    srl_F = [srl_F;mr_F(enl)];
    srr_F = [mr_F(exl);srr_F];
    
    sbl = [sbl;mb(enl)];
    sbr = [mb(exl)*erb;sbr];
    mb(exl) = mb(exl)*(1-erb);
    sbl_C = [sbl_C;mb_C(enl)];
    sbr_C = [mb_C(exl);sbr_C];
    sbl_F = [sbl_F;mb_F(enl)];
    sbr_F = [mb_F(exl);sbr_F];
    
    sll = [sll;ml(enl)];
    slr = [ml(exl)*erl;slr];
    ml(exl) = ml(exl)*(1-erl);
    sll_C = [sll_C;ml_C(enl)];
    slr_C = [ml_C(exl);slr_C];
    sll_F = [sll_F;ml_F(enl)];
    slr_F = [ml_F(exl);slr_F];
    
%% Calculations
    [mt_1,stl_1,str_1] = GUI_Roundabout_Sub_CF(mt,stl,str,mt_C,stl_C,str_C,mt_F,stl_F,str_F,enl,exl,w,v);
    [mr_1,srl_1,srr_1] = GUI_Roundabout_Sub_CF(mr,srl,srr,mr_C,srl_C,srr_C,mr_F,srl_F,srr_F,enl,exl,w,v);
    [mb_1,sbl_1,sbr_1] = GUI_Roundabout_Sub_CF(mb,sbl,sbr,mb_C,sbl_C,sbr_C,mb_F,sbl_F,sbr_F,enl,exl,w,v);
    [ml_1,sll_1,slr_1] = GUI_Roundabout_Sub_CF(ml,sll,slr,ml_C,sll_C,slr_C,ml_F,sll_F,slr_F,enl,exl,w,v);

end

function [m_t1,sl_t1,sr_t1] = GUI_Roundabout_Sub_CF(m,sl,sr,m_C,sl_C,sr_C,m_F,sl_F,sr_F,enl,exl,w,v)

    Sm = min([m_F,m],[],2); % Send and recieve main road top piece
    Rm = max(min([m_F,(w/v)*(m_C-m)],[],2),0);
    
    ym = min([Sm(1:end-1),Rm(2:end)],[],2);
    
    sl(end) = sl(end) + ym(enl-1);
    Ssl = min([sl_F,sl],[],2); % Send and recieve for left and right lanes of the incoming and outgoing side roads
    Rsl = max(min([sl_F,(w/v)*(sl_C-sl)],[],2),0);
    Ssr = min([sr_F,sr],[],2);
    Rsr = max(min([sr_F,(w/v)*(sr_C-sr)],[],2),0);
    
    ysl = min([Ssl(1:end-1),Rsl(2:end)],[],2); % Calculate veicle flow
    ysr = min([Ssr(1:end-1),Rsr(2:end)],[],2);
    
    % Calculate new time step
    m_t1 = m(2:end-1) + ym(1:end-1) - ym(2:end);
    m_t1(exl-1) = m_t1(exl-1) + sr(1) - ysr(1);
    m_t1(enl-1) = m_t1(enl-1) + ysl(end);
    
    sl_t1 = sl(2:end-1) + ysl(1:end-1) - ysl(2:end);
    sr_t1 = sr(2:end-1) + ysr(1:end-1) - ysr(2:end);

end