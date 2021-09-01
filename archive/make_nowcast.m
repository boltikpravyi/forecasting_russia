%%% Nowcasting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script produces a nowcast of VV real GRP growth
% using the estimated parameters from a dynamic factor model.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Clear workspace and set paths.
close all; clear; clc;
addpath('functions');

%% User inputs.
series = 'Real GDP: Russia' ; % Nowcasting real GRP
period = '2016q1'; % Forecasting target quarter
file_xls = 'data/VV/2016-03-01.xls'; % new vintage xls file
file_mat = 'data/VV/mat/2016-03-01.mat'; % new vintage mat file

%% Load model specification and first vintage of data.
% Load model specification structure `Spec`
Spec = load_spec('Spec.xlsx');

%% Load DFM estimation results structure `Res`.
Res = load('ResDFM'); % example_DFM.m used the first vintage of data for estimation
%Res.Res.Mx(1,44) = 1.5/4;
%Res.Res.Mx(1,41) = 1.5/12;

%% Update nowcast and decompose nowcast changes into news.

%%% Nowcast update %%%
vintage_old = '2016-02-01'; datafile_old = fullfile('data', 'VV',[vintage_old '.xls']);
vintage_new = '2016-03-01'; datafile_new = fullfile('data', 'VV',[vintage_new '.xls']);

% Load new vintage in a .mat file
[Z, txt, raw] = xlsread(file_xls);
Time = txt(:,1);
Time(1) = [];
Time = datenum(Time, 'dd.mm.yyyy');
Mnem = txt(1,:);
Mnem(1) = [];
clear x txt raw
save(file_mat, 'Z', 'Time', 'Mnem');

% Load datasets for each vintage
[X_old,~   ] = load_data(datafile_old, Spec);
[X_new,Time] = load_data(datafile_new, Spec);

% check if spec used in estimation is consistent with the current spec
if isequal(Res.Spec, Spec)
    Res = Res.Res;
else
    threshold = 1e-6; % Set to 1e-5 for robust estimates
    Res_RU = dfm(X_new, Spec, threshold);
    save('ResDFM', 'Res', 'Spec');
end 

update_nowcast(X_old, X_new, Time, Spec, Res, series, period, vintage_old, vintage_new);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%