%% Processing ERCOT Grid prices data
% This code processes ERCOT's LMP price historical data files, and provides
% a 24-hour grid price profile based ont the average of every hour of all 
% LMPs in the file 

clc
clear

filename = 'ERCOT_GridPrice0620.csv';% Enter ERCOT LMP .csv file name

Times = readcell(filename,'Range','B2:B362065');
Prices = readmatrix(filename,'Range','D2:D362065');

% Obtaining hourly LMP averages for all locations
flag_init = 1;
Avg_cost = [];
for c = 1:1:length(Prices)-1
    L = Times{c} == Times{c+1};
    if mean(L) ~= 1
        flag_end = c;
        Avg_cost = [Avg_cost;(Prices(flag_init:flag_end))'];
        flag_init = c+1;
    end
    if c+1 == length(Prices)
        Avg_cost = [Avg_cost;(Prices(flag_init:length(Prices)))'];
    end
    
end

Grid = mean(Avg_cost,2);

save('GRID.mat','Grid')% Saving as a .mat file for easy loading into a MATLAB script

