%% ECE 6327 Smart Grid Systems Class Project
% Generating parameters for:
% Gen. units, Main grid, Power storage, Wind turbines, Solar panels,
% Main water system, Water treatment plants and reservoirs, Water
% storage.
clear
clc

% Change parameters as needed and run to update PARAM.mat

% Execute GridData_Processing.m to obtain GRID.mat first
load('GRID.mat')

% Select whether Energy management only, Water managemen only, or Energy-Water management
mode = menu('Select Mode','EM with water from main system','EM only','WM only','Water-Energy');
if mode == 0
    error("No selection made.")
end

% GEN UNIT PARAMETERS
% Lenght of vector denotes number of generators
if mode == 3
    gen_min = 0;%[kW]
    gen_max = 0;%[kW]
    C_G = 0;%[$/kWh]
    NL_G = 0;%[$]
    SU_G = 0;%[$]
    u_G_init = 0;
else
    % Generator parameters data
    gen_min = [5 10 12 15 15 20 30];%[kW]
    gen_max = [95 100 120 130 145 180 200];%[kW]
    C_G = [0.04 0.08 0.14 0.2 0.25 0.3 0.33];%[$/kWh]
    NL_G = [1 1.2 1.6 1.8 2.5 3 3.2];%[$]
    SU_G = [1.8 2.5 3.1 4 4.4 5 5.4];%[$]
    u_G_init = zeros(1,length(C_G));
end
    
% MAIN GRID PARAMETERS
Grid_limit = 300*ones(24,1);%[kW]
C_grid_plus = Grid/100;%[$/kWh]
C_grid_minus = (1-.2)*C_grid_plus;%[$/kWh] ; assuming selling price to be 80% of purchasing price

% POWER STORAGE PARAMETERS
% Lenght of vector denotes number of storage units
if mode == 3
    E_min = 0;%[kWh]
    E_max = 0;%[kWh]
    E_init = 0;% Initial charge [kWh]
    PE_limit = 0;%[kW]
else
    % Energy storage parameters data
    E_min = [20 20 10 15 20 25 30 35 35];%[kWh]
    E_max = [120 115 105 120 130 155 170 200 220];%[kWh]
    E_init = E_min;% Initial charge [kWh]
    PE_limit = [50 45 35 50 65 75 90 95 100];%[kW]
end

% WASTE WATER TREATMENT PARAMETERS
% Lenght of vector denotes number of waste water treating units
W_WW = [370];%[gal/kWh]
NL_WW = [1];%[$]
WL_WW_limit = [25000];%[gal]
if mode == 3 || mode == 4
    W_WW_min = [20];%[gal/h]
    W_WW_max = [4000];%[gal/h]
    WL_WW_init = 0.1*WL_WW_limit;%[gal]
    reclaimed_w = 50;% percent of water demand reclaimed as waste water [%]
    S_WW = 35;% Surface area of wastewater reservoir[m^2]
    WR_init = zeros(24,1);
else
    % Energy management only
    W_WW_min = 0;%[gal/h]
    W_WW_max = 0;%[gal/h]
    WL_WW_init = 0*WL_WW_limit;
    reclaimed_w = 0;% percent of water demand reclaimed as waste water [%]
    S_WW = 0;% Surface area of wastewater reservoir[m^2]
    WR_init = 0;
end

% MAIN WATER SYSTEM PARAMETERS
if mode == 2
    Wmain_limit = 0;%[gal/h]
    C_Wmain_plus = 0;%[$/gal]
    C_Wmain_minus = 0;%[$/gal]
else
    Wmain_limit = 3000*ones(24,1);%[gal/h]
%     Wmain_limit(Fault) = 0;
    C_Wmain_plus = 0.01;%[$/gal]
    C_Wmain_minus = (1-.2)*C_Wmain_plus;%[$/gal]
end

% WATER STORAGE PARAMETERS
% Lenght of vector denotes number of water storage units
if mode == 3 || mode == 4
    WS_limit = [350 550 700 850];%[gal/h]
    WL_ST_limit = [10000 16000 21000 30000];%[gal]
    WL_ST_init = zeros(size(WS_limit));%[gal]
else
    WS_limit = 0;%[gal/h]
    WL_ST_limit = 0;%[gal]
    WL_ST_init = 0;%[gal]
end

if mode == 3
    % SOLAR POWER
    % Lenght of vector denotes number of solar panels
    S_SP = 0;%[m^2]
    effi_SP = 0;%[-]
    tilt_a = 0;%[°]
    collect_a = 0;%[°]
    % WIND POWER
    % Lenght of vector denotes number of wind turbines
    Prated_WT = 0;%[kW]
    CutIn = 0;%[m/s]
    v_rated = 0;%[m/s]
    CutOut = 0;%[m/s]
else
    % SOLAR POWER
    % Lenght of vector denotes number of solar panels
    S_SP = [(4*4)*ones(1,25) (4*4)*ones(1,25) (3*3)*ones(1,20) (6*5)*ones(1,4)];%[m^2]
    effi_SP = [0.35*ones(1,25) 0.35*ones(1,25) 0.3*ones(1,20) 0.4*ones(1,4)];%[-]
    tilt_a = [25*ones(1,25) 20*ones(1,25) 40*ones(1,20) 35*ones(1,4)];%[°]
    collect_a =[45*ones(1,25) -45*ones(1,25) 0*ones(1,20) 0*ones(1,4)];%[°]
    % WIND POWER
    % Lenght of vector denotes number of wind turbines
    Prated_WT = [100 150];%[kW]
    CutIn = [2.5 3.2];%[m/s]
    v_rated = [10 14];%[m/s]
    CutOut = [13.2  20];%[m/s]
end

% Total number of loads
Rload_num = 60;% residential units
Cload_num = 2;% commercial units

% Load shedding cost
C_shedW = 100000;%[$/gal]
C_shedP = 100000;%[$/kWh]

save('PARAM.mat')
clear
clc

