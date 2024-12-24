%% Water-Energy Management Co-Optimization for a Community Microgrid
%  Jesus Silva Rodriguez | UH Cullen College of Engineering
% ***Ensure Gurobi is enabled. Enter in command window: cvx_solver gurobi_2***
clear
clc
close all
%% Data sets
% LOADING DATA
SOLAR = readmatrix('SOLAR.csv','Range','C2163:M2882');
WIND = readmatrix('WIND_TEMP.xlsx','Range','D11:D346');
RESIDENT_LOAD = readmatrix('RESIDENT_LOAD.csv','Range','B2162:M2881');
COMM_LOAD = readmatrix('COMM_LOAD.csv','Range','B2162:M2881');
dates_data = readcell('RAIN2.xlsx','Range','F2:F79');
RAIN = readmatrix('RAIN2.xlsx','Range','G2:G79');
RESIDENT_WATER = readmatrix('WATER_D.xlsx','Range','C3:C26');
COMM_WATER = readmatrix('WATER_D.xlsx','Range','D3:D26');

% PROCESSING DATA
% Solar irradiance
GHI = [];%[W/m^2]
DNI = [];%[W/m^2]
DHI = [];%[W/m^2]
for c = 1:24:length(SOLAR)
    GHI = [GHI SOLAR(c:c+23,3)];
    DNI = [DNI SOLAR(c:c+23,6)];
    DHI = [DHI SOLAR(c:c+23,9)];
end
GHI_av = mean(GHI,2);% average GHI at every time interval [W/m^2]
DNI_av = mean(DNI,2);% average DNI at every time interval [W/m^2]
DHI_av = mean(DHI,2);% average DHI at every time interval [W/m^2]

% Wind speeds
WS = [];%[km/h]
for c = 1:24:length(WIND)
    WS = [WS WIND(c:c+23)];
end
v_av = (1000/3600)*mean(WS,2);% average wind speeds at every time interval [m/s]

% Rain
Dates = [];
for d = 1:1:30
    for h = 0:1:23
    D = sprintf("%d-Apr-2013 %d:00:00",d,h);
    Dates = [Dates;D];
    end
end
Dates = datetime(Dates);
for i = 1:1:24*30 - 1
    if Dates(i) ~= dates_data{i}
        dates_data(i+1:end+1,:) = dates_data(i:end,:);
        RAIN(i+1:end+1,:) = RAIN(i:end,:);
        dates_data{i,1} = Dates(i);
        RAIN(i,1) = 0;
    elseif dates_data{i} == dates_data{end}
            dates_data{i+1,1} = Dates(i+1);
            RAIN(i+1,1) = 0; 
    end
end
Rain = [];%[mm/h]
for c = 1:24:length(RAIN)
    Rain = [Rain RAIN(c:c+23)];
end
Rh_av = (1/1000)*mean(Rain,2);% average precipitation at every time interval [m/h]

% Residential and commercial power loads
RESIDENT_LOAD = sum(RESIDENT_LOAD,2);
COMM_LOAD = sum(COMM_LOAD,2);
Res_L = [];% [kW]
Comm_L = [];% [kW]
for c = 1:24:length(COMM_LOAD)
    Res_L = [Res_L RESIDENT_LOAD(c:c+23)];
    Comm_L = [Comm_L COMM_LOAD(c:c+23)];
end
LR_av = mean(Res_L,2);% average residential load at every time interval [kW]
LC_av = mean(Comm_L,2);% average commercial load at every time interval [kW]

% Residential and commercial water loads
DR_av = RESIDENT_WATER;% average residential water demand at every time interval [m^3/h]
DC_av = COMM_WATER;% average commercial water demand at every time interval [m^3/h]

% LOADING PARAMETERS
load('PARAM.mat')

%% Renewable Energy
% SOLAR POWER
Lat = readmatrix('SOLAR.csv','Range','E1:E1');%[°]
n = (90+120)/2;% A day in the middle of April
delta = 23.45*sind((360/365)*(n-81));%[°]
Hour_a = 15*(12-[1:24])';%[°]
Alt_a = rad2deg(asin(cosd(Lat)*cosd(delta)*cosd(Hour_a) + sind(Lat)*sind(delta)));%[°]
Solar_a = rad2deg(asin((cosd(delta)*sind(Hour_a))./cosd(Alt_a)));%[°]
incidence = [];%[rad]
% Incidence angle- Rows: time intervals ; Columns: solar panels
for t = 1:1:24
    incidence = [incidence;acos(cosd(Alt_a(t))*cosd(Solar_a(t) - collect_a).*sind(tilt_a) + sind(Alt_a(t))*cosd(tilt_a))];
end
incidence = rad2deg(incidence);%[°]
reflectance = 0.2;% Ground light reflectance [-]

% Calculating irradiances
I_BCav = [];%[W/m^2]
I_DCav = [];%[W/m^2]
I_RCav = [];%[W/m^2]
for t = 1:1:24
    I_BCav = [I_BCav;DNI_av(t)*cosd(incidence(t,:))];
    I_DCav = [I_DCav;DHI_av(t)*((1+cosd(tilt_a))/2)];
    I_RCav = [I_RCav;GHI_av(t)*reflectance*((1-cosd(tilt_a))/2)];
end
I_av = I_BCav + I_DCav + I_RCav;%[W/m^2]

% Power
P_SP = [];%[W]
for t = 1:1:24
    P_SP = [P_SP;effi_SP.*S_SP.*I_av(t,:)];
end
P_SP(P_SP<0) = 0;
P_SP = P_SP/1000;%[kW]

% WIND POWER
P_WT = [];%[kW]
for t = 1:1:24
    p_wt = Prated_WT.*((v_av(t)-CutIn)./(v_rated-CutIn));
    for s = 1:1:length(Prated_WT)
        if v_av(t)<CutIn(s) || v_av(t)>CutOut(s)
            p_wt(s) = 0;
        elseif v_av(t)>v_rated(s) && v_av(t)<=CutOut(s)
            p_wt(s) = Prated_WT;
        end
    end
    P_WT = [P_WT;p_wt];
end

%% Water Supply
% RAIN WATER
% Rain water obtained per hour
RR = [];%[m^3/h]
for t = 1:1:24
    RR = [RR;Rh_av(t)*S_WW];
end

% WASTE WATER
% Waste water obtained per hour
WR = zeros(1,length(W_WW));%[m^3/h]
for t = 2:1:24
    WR = [WR;((reclaimed_w/100)*((Rload_num*DR_av(t-1) + Cload_num*DC_av(t-1))/length(W_WW)))];
end
WR = (WR + RR)*264.172;% Combining rain and waste water into same untreated reservoir [gal/h]

DR_av = DR_av*264.172;% [gal/h]
DC_av = DC_av*264.172;% [gal/h]
%% Optimization

% EM with water from main system Case
if mode == 1
    
cvx_begin
% VARIABLES
    % Binary Unit Commitment variables
    variable u_G(length(NL_G),24) binary
    variable v_G(length(SU_G),24) binary
    variable p_plus(1,24) binary
    variable p_minus(1,24) binary
    variable e_c(length(PE_limit),24) binary
    variable e_d(length(PE_limit),24) binary
    variable a_plus(1,24) binary
    % Power variables
    variables P_G(length(gen_min),24) P_grid_plus(1,24) P_grid_minus(1,24)
    variables PE_d(length(PE_limit),24) PE_c(length(PE_limit),24) E(length(PE_limit),24)
    % Water variables
    variable Wmain_plus(1,24)
    
% OBJECTIVE FUNCTION
    % Unit commitment
    f_U = sum(NL_G*u_G + SU_G*v_G);
    f_E = sum(C_G*P_G + (C_grid_plus').*P_grid_plus - (C_grid_minus').*P_grid_minus);
    f_W = sum((C_Wmain_plus').*Wmain_plus);
    minimize f_U + f_E + f_W
    
% CONSTRAINTS
    subject to
        % Generator Constraints
        for t = 1:1:24
            (gen_min').*u_G(:,t) <= P_G(:,t) <= (gen_max').*u_G(:,t)
        end
        v_G(:,1) >= u_G(:,1)
        for t = 2:1:24
            v_G(:,t) >= u_G(:,t) - u_G(:,t-1)
        end
        % Main grid constraints
        0 <= P_grid_plus <= (Grid_limit').*p_plus
        0 <= P_grid_minus <= (Grid_limit').*p_minus
        % Power storage constraints
        for t = 1:1:24
            0 <= PE_d(:,t) <= (PE_limit').*e_d(:,t)
            0 <= PE_c(:,t) <= (PE_limit').*e_c(:,t)
            (E_min') <= E(:,t) <= (E_max')
        end
        e_c + e_d <= 1
        E(:,1) == E_init' + (PE_c(:,1) - PE_d(:,1))
        for t = 2:1:24
            E(:,t) == E(:,t-1) + (PE_c(:,t) - PE_d(:,t))
        end
        E(:,end) == E_init'
        % Power Balance
        sum(P_G) + P_grid_plus - P_grid_minus + sum(PE_d - PE_c) == Rload_num*(LR_av') + Cload_num*(LC_av') - sum(P_SP') - sum(P_WT')
        % Main Water System
        0 <= Wmain_plus <= (Wmain_limit').*a_plus
        % Water Balance
        Wmain_plus == Rload_num*(DR_av') + Cload_num*(DC_av')
        
cvx_end
clc

% EM Only
elseif mode == 2
    
cvx_begin
% VARIABLES
    % Binary Unit Commitment variables
    variable u_G(length(NL_G),24) binary
    variable v_G(length(SU_G),24) binary
    variable p_plus(1,24) binary
    variable p_minus(1,24) binary
    variable e_c(length(PE_limit),24) binary
    variable e_d(length(PE_limit),24) binary
    % Power variables
    variables P_G(length(gen_min),24) P_grid_plus(1,24) P_grid_minus(1,24)
    variables PE_d(length(PE_limit),24) PE_c(length(PE_limit),24) E(length(PE_limit),24)
    
% OBJECTIVE FUNCTION
    % Unit commitment
    f_U = sum(NL_G*u_G + SU_G*v_G);
    % Economic dispatch
    f_E = sum(C_G*P_G + (C_grid_plus').*P_grid_plus - (C_grid_minus').*P_grid_minus);
    minimize f_U + f_E
    
% CONSTRAINTS
    subject to
        % Generator Constraints
        for t = 1:1:24
            (gen_min').*u_G(:,t) <= P_G(:,t) <= (gen_max').*u_G(:,t)
        end
        v_G(:,1) >= u_G(:,1)
        for t = 2:1:24
            v_G(:,t) >= u_G(:,t) - u_G(:,t-1)
        end
        % Main grid constraints
        0 <= P_grid_plus <= (Grid_limit').*p_plus
        0 <= P_grid_minus <= (Grid_limit').*p_minus
        % Power storage constraints
        for t = 1:1:24
            0 <= PE_d(:,t) <= (PE_limit').*e_d(:,t)
            0 <= PE_c(:,t) <= (PE_limit').*e_c(:,t)
            (E_min') <= E(:,t) <= (E_max')
        end
        e_c + e_d <= 1
        E(:,1) == E_init' + (PE_c(:,1) - PE_d(:,1))
        for t = 2:1:24
            E(:,t) == E(:,t-1) + (PE_c(:,t) - PE_d(:,t))
        end
        E(:,end) == E_init'
        % Power Balance
        sum(P_G) + P_grid_plus - P_grid_minus + sum(PE_d - PE_c) == Rload_num*(LR_av') + Cload_num*(LC_av') - sum(P_SP') - sum(P_WT')

cvx_end
clc

% WM Only
elseif mode == 3

cvx_begin
%VARIABLES
    % Binary Unit Commitment variables
    variable u_WW(length(NL_WW),24) binary
    variable a_plus(1,24) binary
    variable a_minus(1,24) binary
    variable r_c(length(WS_limit),24) binary
    variable r_d(length(WS_limit),24) binary
    % Water variables
    variables L_WW(length(W_WW),24) WL_WW(length(W_WW),24) WL_ST(length(WS_limit),24)
    variables WS_c(length(WS_limit),24) WS_d(length(WS_limit),24) Wmain_plus(1,24) Wmain_minus(1,24)
    
% OBJECTIVE FUNCTION
    % Unit commitment
    f_U = sum(NL_WW*u_WW);
    % Economic dispatch
    f_W = sum((C_grid_plus').*L_WW + (C_Wmain_plus').*Wmain_plus - (C_Wmain_minus').*Wmain_minus);
    
    minimize f_U + f_W
    
% CONSTRAINTS
    subject to
        % Water treatment constraints
        for t = 1:1:24
            (W_WW_min').*u_WW(:,t) <= (W_WW').*L_WW(:,t) <= (W_WW_max').*u_WW(:,t)
            0 <= WL_WW(:,t) <= (WL_WW_limit')
        end
        WL_WW(:,1) == (WL_WW_init') + ((WR(1,:)') - (W_WW').*L_WW(:,1))
        for t = 2:1:24
            WL_WW(:,t) == WL_WW(:,t-1) + ((WR(t,:)') - (W_WW').*L_WW(:,t))
        end
        % Main water system constraints
        0 <= Wmain_plus <= (Wmain_limit').*a_plus
        0 <= Wmain_minus <= (Wmain_limit').*a_minus
        a_plus + a_minus <= 1
        % Water storage constraints
        for t = 1:1:24
            0 <= WS_c(:,t) <= (WS_limit').*r_c(:,t)
            0 <= WS_d(:,t) <= (WS_limit').*r_d(:,t)
        end
        r_c + r_d <= 1
        WL_ST(:,1) == (WL_ST_init') + (WS_c(:,1) - WS_d(:,1))
        for t = 2:1:24
            WL_ST(:,t) == WL_ST(:,t-1) + (WS_c(:,t) - WS_d(:,t))
            0 <= WL_ST(:,t) <= (WL_ST_limit')
        end
        WL_ST(:,end) == WL_ST_init'
        % Water balance
        W_WW*L_WW + Wmain_plus - Wmain_minus + sum(WS_d - WS_c) == Rload_num*(DR_av') + Cload_num*(DC_av')
        % Main grid
        0 <= L_WW <= Grid_limit'
        
cvx_end
clc

% Water-Energy
elseif mode == 4

cvx_begin
% VARIABLES
    % Binary Unit Commitment variables
    variable u_G(length(NL_G),24) binary
    variable v_G(length(SU_G),24) binary
    variable u_WW(length(NL_WW),24) binary
    variable p_plus(1,24) binary
    variable p_minus(1,24) binary
    variable e_c(length(PE_limit),24) binary
    variable e_d(length(PE_limit),24) binary
    variable a_plus(1,24) binary
    variable a_minus(1,24) binary
    variable r_c(length(WS_limit),24) binary
    variable r_d(length(WS_limit),24) binary
    % Power variables
    variables P_G(length(gen_min),24) P_grid_plus(1,24) P_grid_minus(1,24)
    variables PE_d(length(PE_limit),24) PE_c(length(PE_limit),24) E(length(PE_limit),24)
    % Water variables
    variables L_WW(length(W_WW),24) WL_WW(length(W_WW),24) WL_ST(length(WS_limit),24)
    variables WS_c(length(WS_limit),24) WS_d(length(WS_limit),24) Wmain_plus(1,24) Wmain_minus(1,24)
    
% OBJECTIVE FUNCTION
    % Unit commitment
    f_U = sum(NL_G*u_G + SU_G*v_G + NL_WW*u_WW);
    % Economic dispatch
    f_E = sum(C_G*P_G + (C_grid_plus').*P_grid_plus - (C_grid_minus').*P_grid_minus);
    f_W = sum((C_Wmain_plus').*Wmain_plus - (C_Wmain_minus').*Wmain_minus);
    
    minimize f_U + f_E + f_W
    
% CONSTRAINTS
    subject to
        % Generator constraints
        for t = 1:1:24
            (gen_min').*u_G(:,t) <= P_G(:,t) <= (gen_max').*u_G(:,t)
        end
        v_G(:,1) >= u_G(:,1)
        for t = 2:1:24
            v_G(:,t) >= u_G(:,t) - u_G(:,t-1)
        end
        % Main grid constraints
        0 <= P_grid_plus <= (Grid_limit').*p_plus
        0 <= P_grid_minus <= (Grid_limit').*p_minus
        p_plus + p_minus <= 1
        % Power storage constraints
        for t = 1:1:24
            0 <= PE_d(:,t) <= (PE_limit').*e_d(:,t)
            0 <= PE_c(:,t) <= (PE_limit').*e_c(:,t)
            (E_min') <= E(:,t) <= (E_max')
        end
        e_c + e_d <= 1
        E(:,1) == E_init' + (PE_c(:,1) - PE_d(:,1))
        for t = 2:1:24
            E(:,t) == E(:,t-1) + (PE_c(:,t) - PE_d(:,t))
        end
        E(:,end) == E_init'
        % Power balance
        sum(P_G) + P_grid_plus - P_grid_minus + sum(PE_d - PE_c) + sum(P_SP') + sum(P_WT') == Rload_num*(LR_av') + Cload_num*(LC_av') + L_WW
        % Water treatment constraints
        for t = 1:1:24
            (W_WW_min').*u_WW(:,t) <= (W_WW').*L_WW(:,t) <= (W_WW_max').*u_WW(:,t)
            0 <= WL_WW(:,t) <= (WL_WW_limit')
        end
        WL_WW(:,1) == (WL_WW_init') + ((WR(1,:)') - (W_WW').*L_WW(:,1))
        for t = 2:1:24
            WL_WW(:,t) == WL_WW(:,t-1) + ((WR(t,:)') - (W_WW').*L_WW(:,t))
        end
        % Main water system constraints
        0 <= Wmain_plus <= (Wmain_limit').*a_plus
        0 <= Wmain_minus <= (Wmain_limit').*a_minus
        a_plus + a_minus <= 1
        % Water storage constraints
        for t = 1:1:24
            0 <= WS_c(:,t) <= (WS_limit').*r_c(:,t)
            0 <= WS_d(:,t) <= (WS_limit').*r_d(:,t)
        end
        r_c + r_d <= 1
        WL_ST(:,1) == (WL_ST_init') + (WS_c(:,1) - WS_d(:,1))
        for t = 2:1:24
            WL_ST(:,t) == WL_ST(:,t-1) + (WS_c(:,t) - WS_d(:,t))
            0 <= WL_ST(:,t) <= (WL_ST_limit')
        end
        WL_ST(:,end) == WL_ST_init'
        % Water balance
        W_WW*L_WW + Wmain_plus - Wmain_minus + sum(WS_d - WS_c) == Rload_num*(DR_av') + Cload_num*(DC_av')
        
cvx_end
clc

end
  
%% Results
T = 1:24;
fprintf('\tRESULTS\n\nObj. Value: $%0.2f\n',cvx_optval)

% ELECTRICITY RESULTS
if mode == 1 || mode == 2
    fprintf('Energy Cost: $%0.2f\n',f_E+f_U)
    
    L_net = Rload_num*(LR_av') + Cload_num*(LC_av') - sum(P_SP') - sum(P_WT');
    PE_netC = sum(PE_c) - sum(PE_d);
    PE_netC(PE_netC<0) = 0;
    PE_netD = sum(PE_d) - sum(PE_c);
    PE_netD(PE_netD<0) = 0;
    

    figure(1)
    stairs(T,L_net,'color',[0 0.5 0],'LineWidth',0.8)
    hold on
    stairs(T,sum(P_G),'color',[1 0.647 0],'LineWidth',0.8)
    stairs(T,P_grid_plus,'r','LineWidth',0.8)
    stairs(T,P_grid_minus,'r--','LineWidth',0.8)
    stairs(T,PE_netC,'b--','LineWidth',0.8)
    stairs(T,PE_netD,'b','LineWidth',0.8)
    ax = gca;
    ax.FontSize = 13;
    axis([T(1) T(end) 0 (1+.1)*max(L_net)])
    ylabel('Power [kW]')
    xlabel('Time of Day [hr]')
    legend('Net Load','Gen.','From Grid','To Grid','ES Net Charge','ES Net Discharge','Location', 'NorthWest','FontSize',10)

    hold off
    figure(2)
    yyaxis left
    stairs(T,sum(E),'b','LineWidth',0.8)
    hold on
    yyaxis right
    stairs(T,C_grid_plus','-','LineWidth',0.8)
    ax = gca;
    ax.FontSize = 13;
    axis([T(1) T(end) 0 inf])
    yyaxis left
    ylabel('Energy Stored [kWh]')
    yyaxis right
    ylabel('Electricity Price [$/kWh]')
    xlabel('Time of Day [hr]')
    legend('Energy Storage','Main Grid Price','location','Northeast','FontSize',10)
    hold off

elseif mode == 3
    
    figure(1)
    stairs(T,L_WW,'r','LineWidth',0.8)
    axis([T(1) T(end) 0 (1+.1)*max(L_WW)])
    ax = gca;
    ax.FontSize = 13;
    ylabel('Power [kW]')
    xlabel('Time of Day [hr]')
    legend('From Grid','location','Northwest','FontSize',10)
    f = 2;

elseif mode == 4
    fprintf('Energy Cost: $%0.2f\n',f_E+sum(NL_G*u_G + SU_G*v_G))
    
    L_net = Rload_num*(LR_av') + Cload_num*(LC_av') + L_WW - sum(P_SP') - sum(P_WT');
    PE_netC = sum(PE_c) - sum(PE_d);
    PE_netC(PE_netC<0) = 0;
    PE_netD = sum(PE_d) - sum(PE_c);
    PE_netD(PE_netD<0) = 0;
    
    figure(1)
    stairs(T,L_net,'color',[0 0.5 0],'LineWidth',0.8)
    hold on
    stairs(T,sum(P_G),'color',[1 0.647 0],'LineWidth',0.8)
    stairs(T,P_grid_plus,'r','LineWidth',0.8)
    stairs(T,P_grid_minus,'r--','LineWidth',0.8)
    stairs(T,PE_netC,'b--','LineWidth',0.8)
    stairs(T,PE_netD,'b','LineWidth',0.8)
    ax = gca;
    ax.FontSize = 13;
    axis([T(1) T(end) 0 (1+.1)*max(L_net)])
    ylabel('Power [kW]')
    xlabel('Time of Day [hr]')
    legend('Net Load','Gen.','From Grid','To Grid','ES Net Charge','ES Net Discharge','Location', 'NorthWest','FontSize',10)

    hold off
    f = 3;
    
    figure(2)
    yyaxis left
    stairs(T,sum(E),'b','LineWidth',0.8)
    hold on
    yyaxis right
    stairs(T,C_grid_plus','-','LineWidth',0.8)
    ax = gca;
    ax.FontSize = 13;
    axis([T(1) T(end) 0 inf])
    yyaxis left
    ylabel('Energy Stored [kWh]')
    yyaxis right
    ylabel('Electricity Price [$/kWh]')
    xlabel('Time of Day [hr]')
    legend('Energy Storage','Main Grid Price','location','Northeast','FontSize',10)
    hold off
    
end

% WATER RESULTS
TWD = Rload_num*(DR_av') + Cload_num*(DC_av');

if mode == 1
    fprintf('Water Cost: $%0.2f\n',f_W)
    
    figure(3)
    stairs(T,TWD,'color',[0 0.5 0],'LineWidth',0.8)
    ax = gca;
    ax.FontSize = 13;
    ylabel('Water Flow [gal/h]')
    xlabel('Time of Day [hr]')
    legend('Water Demand','location','Northwest','FontSize',10)
    axis([T(1) T(end) 0 (1+.1)*max(TWD)])

elseif mode == 3 || mode == 4
    figure(f)
    fprintf('Water Cost: $%0.2f\n',f_W + sum(NL_WW*u_WW))
    
    WS_netC = sum(WS_c) - sum(WS_d);
    WS_netC(WS_netC<0) = 0;
    WS_netD = sum(WS_d) - sum(WS_c);
    WS_netD(WS_netD<0) = 0;
    
    stairs(T,TWD,'color',[0 0.5 0],'LineWidth',0.8)
    hold on
    stairs(T,W_WW*L_WW,'color',[1 0.647 0],'LineWidth',0.8)
    stairs(T,Wmain_plus,'r','LineWidth',0.8)
    stairs(T,Wmain_minus,'r--','LineWidth',0.8)
    stairs(T,WS_netC,'b--','LineWidth',0.8)
    stairs(T,WS_netD,'b','LineWidth',0.8)
    ax = gca;
    ax.FontSize = 13;
    ylabel('Water Flow [gal/h]')
    xlabel('Time of Day [hr]')
    legend('Demand','Treated','From Main','To Main','Net Stored','Net Released','location','Northeast','FontSize',10)
    if max(W_WW*L_WW) > max(TWD)
        axis([T(1) T(end) 0 (1+.1)*max(W_WW*L_WW)])
    else
        axis([T(1) T(end) 0 (1+.1)*max(TWD)])
    end

end
