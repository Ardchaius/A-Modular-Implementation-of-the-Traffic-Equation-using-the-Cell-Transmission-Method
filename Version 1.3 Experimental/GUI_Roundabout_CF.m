% Improved Roundabout Segment function

function [mt_1,mr_1,mb_1,ml_1,stl_1,str_1,srl_1,srr_1,sbl_1,sbr_1,sll_1,slr_1]...
    = GUI_Roundabout_CF(mt,mr,mb,ml,stl,str,srl,srr,sbl,sbr,sll,slr,...
                    mt_C,mr_C,mb_C,ml_C,stl_C,str_C,srl_C,srr_C,sbl_C,sbr_C,sll_C,slr_C,...
                    mt_F,mr_F,mb_F,ml_F,stl_F,str_F,srl_F,srr_F,sbl_F,sbr_F,sll_F,slr_F,length)

%% Constants

    exl = length/2 - 0.5; % Location of exit lane in direction of travel adjusted for input and output cells
    enl = length/2 + 1.5; % Location of entrance lane in direction of travel adjusted for input and output cells
    ert = 0.25; % Exit rate of top piece
    err = 0.25; % Exit rate of right piece
    erb = 0.25; % Exit rate of bottom piece
    erl = 0.25; % Exit rate of left piece
    
%% Initializations
    
    % Format data to correct format for use in calculations
    splitVeh = SplitFunc(mt(exl),ert);
    stl = [stl;mt(enl)];
    str = [splitVeh(1);str];
    mt(exl) = splitVeh(2);
    stl_C = [stl_C;mt_C(enl)];
    str_C = [mt_C(exl);str_C];
    stl_F = [stl_F;mt_F(enl)];
    str_F = [mt_F(exl);str_F];
    
    splitVeh = SplitFunc(mr(exl),err);
    srl = [srl;mr(enl)];
    srr = [splitVeh(1);srr];
    mr(exl) = splitVeh(2);
    srl_C = [srl_C;mr_C(enl)];
    srr_C = [mr_C(exl);srr_C];
    srl_F = [srl_F;mr_F(enl)];
    srr_F = [mr_F(exl);srr_F];
    
    splitVeh = SplitFunc(mb(exl),erb);
    sbl = [sbl;mb(enl)];
    sbr = [splitVeh(1);sbr];
    mb(exl) = splitVeh(2);
    sbl_C = [sbl_C;mb_C(enl)];
    sbr_C = [mb_C(exl);sbr_C];
    sbl_F = [sbl_F;mb_F(enl)];
    sbr_F = [mb_F(exl);sbr_F];
    
    splitVeh = SplitFunc(ml(exl),erl);
    sll = [sll;ml(enl)];
    slr = [splitVeh(1);slr];
    ml(exl) = splitVeh(2);
    sll_C = [sll_C;ml_C(enl)];
    slr_C = [ml_C(exl);slr_C];
    sll_F = [sll_F;ml_F(enl)];
    slr_F = [ml_F(exl);slr_F];
    
%% Calculations
    [mt_1,stl_1,str_1] = GUI_Roundabout_Sub_CF(mt,stl,str,mt_C,stl_C,str_C,mt_F,stl_F,str_F,enl,exl);
    [mr_1,srl_1,srr_1] = GUI_Roundabout_Sub_CF(mr,srl,srr,mr_C,srl_C,srr_C,mr_F,srl_F,srr_F,enl,exl);
    [mb_1,sbl_1,sbr_1] = GUI_Roundabout_Sub_CF(mb,sbl,sbr,mb_C,sbl_C,sbr_C,mb_F,sbl_F,sbr_F,enl,exl);
    [ml_1,sll_1,slr_1] = GUI_Roundabout_Sub_CF(ml,sll,slr,ml_C,sll_C,slr_C,ml_F,sll_F,slr_F,enl,exl);

end

function [m_t1,sl_t1,sr_t1] = GUI_Roundabout_Sub_CF(m,sl,sr,m_C,sl_C,sr_C,m_F,sl_F,sr_F,enl,exl)

    Sm = min([m_F,m],[],2); % Send and recieve main road top piece
    Rm = max(min([m_F,(m_C-m)],[],2),0);
    
    ym = min([Sm(1:end-1),Rm(2:end)],[],2);
    
    sl(end) = sl(end) + ym(enl-1);
    Ssl = min([sl_F,sl],[],2); % Send and recieve for left and right lanes of the incoming and outgoing side roads
    Rsl = max(min([sl_F,(sl_C-sl)],[],2),0);
    Ssr = min([sr_F,sr],[],2);
    Rsr = max(min([sr_F,(sr_C-sr)],[],2),0);
    
    ysl = min([Ssl(1:end-1),Rsl(2:end)],[],2); % Calculate veicle flow
    ysr = min([Ssr(1:end-1),Rsr(2:end)],[],2);
    
    % Calculate new time step
    m_t1 = m(2:end-1) + ym(1:end-1) - ym(2:end);
    m_t1(exl-1) = m_t1(exl-1) + sr(1) - ysr(1);
    m_t1(enl-1) = m_t1(enl-1) + ysl(end);
    
    sl_t1 = sl(2:end-1) + ysl(1:end-1) - ysl(2:end);
    sr_t1 = sr(2:end-1) + ysr(1:end-1) - ysr(2:end);

end