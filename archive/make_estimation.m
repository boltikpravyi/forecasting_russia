%%% Dynamic factor model (DFM) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script estimates a dynamic factor model (DFM) using a panel of
% monthly and quarterly series.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Clear workspace and set paths.
close all; clear; clc;
addpath('functions');

%% User inputs.
vintage = '2021-08-19'; % vintage dataset to use for estimation
country = 'VV';         % Russian macroeconomic data
sample_start = datenum('2002-01-01', 'yyyy-mm-dd'); % estimation sample

%% Load model specification and dataset.
% Load model specification structure `Spec`
Spec = load_spec('Spec.xlsx');
% Parse `Spec`
SeriesID = Spec.SeriesID; SeriesName = Spec.SeriesName; Units = Spec.Units; UnitsTransformed = Spec.UnitsTransformed;

% Load data
datafile = fullfile('data', country, [vintage '.xls']);
[X, Time, Z] = load_data(datafile, Spec, sample_start);
summarize(X, Time, Spec, vintage); % summarize data

%% Run dynamic factor model (DFM) and save estimation output as 'ResDFM'.
threshold = 1e-6; % Set to 1e-5 for more robust estimates

Res = dfm(X, Spec, threshold);
save('ResDFM', 'Res', 'Spec');

%% Plot projection of common factor.
idxSeries = strcmp('Real GDP: Russia', SeriesID); t_obs = ~isnan(X(:,idxSeries));
figure('Name', 'Common Factor Projection');
CommonFactor = Res.C(idxSeries, 1:5) * Res.Z(:, 1:5)' * Res.Wx(idxSeries) + Res.Mx(idxSeries);
plot(Time, CommonFactor, 'k'); hold on;
plot(Time(t_obs), X(t_obs, idxSeries), 'b'); box on;
title(SeriesName{idxSeries}); xlim(Time([1 end])); datetick('x', 'yyyy', 'keeplimits');
ylabel({Units{idxSeries}, UnitsTransformed{idxSeries}});
legend('common component', 'data'); legend boxoff;

GlobalFactor = Res.C(idxSeries, 1) * Res.Z(:, 1)' * Res.Wx(idxSeries);
RealFactor = Res.C(idxSeries, 11) * Res.Z(:, 11)' * Res.Wx(idxSeries);
FinancialFactor = Res.C(idxSeries, 16) * Res.Z(:, 16)' * Res.Wx(idxSeries);
RegionalFactor = Res.C(idxSeries, 21) * Res.Z(:, 21)' * Res.Wx(idxSeries);
plot(Time, GlobalFactor + RealFactor + FinancialFactor + RegionalFactor + Res.Mx(idxSeries), 'k'); hold on;
plot(Time(t_obs), X(t_obs, idxSeries), 'b'); box on;

idxSeries = strcmp('GRP', SeriesID); t_obs = ~isnan(X(:,idxSeries));
GlobalFactor = Res.C(idxSeries, 1:5) * Res.Z(:, 1:5)' * Res.Wx(idxSeries);
RealFactor = Res.C(idxSeries, 11:15) * Res.Z(:, 11:15)' * Res.Wx(idxSeries);
RegionalFactor = Res.C(idxSeries, 21:25) * Res.Z(:, 21:25)' * Res.Wx(idxSeries);
plot(Time, GlobalFactor, 'k'); hold on;
plot(Time(t_obs), X(t_obs,idxSeries), 'b'); box on;
title(SeriesName{idxSeries}); xlim(Time([1 end])); datetick('x', 'yyyy', 'keeplimits');
ylabel({Units{idxSeries}, UnitsTransformed{idxSeries}});
legend('common component', 'data'); legend boxoff;

%% Idyosyncratic components
y = Res.Z';
y = y(:,2:end);
x = [Res.Z(:,2:5) lagmatrix(Res.Z(:,5),1) Res.Z(:,7:10) lagmatrix(Res.Z(:,10),1) Res.Z(:,12:15) ...
    lagmatrix(Res.Z(:,15),1) Res.Z(:,17:20) lagmatrix(Res.Z(:,20),1) Res.Z(:,22:25) lagmatrix(Res.Z(:,25),1) lagmatrix(Res.Z(:,26:68), 1) ...
    Res.Z(:,70:73) lagmatrix(Res.Z(:,73),1)];
x = x(2:end,:);
x = x';
A = Res.A;
e = (y - A * x)';
