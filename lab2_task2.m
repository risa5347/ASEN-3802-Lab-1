clc; clear; close all;

% File names for each loading case
case1_file = 'case_1_data.txt';
case2_file = 'case_2_data.txt';
case3_file = 'case_3_data.txt';

% Load raw experimental data
C1 = readmatrix(case1_file);
C2 = readmatrix(case2_file);
C3 = readmatrix(case3_file);

% Map data columns to physical measurements
getCols = @(C) struct( ...
    'Load', C(:,1), ...
    'F0',   C(:,2), ...
    'F1',   C(:,3), ...
    'F2',   C(:,4), ...
    'F3D',  C(:,5), ...
    'LVDT', C(:,6));

S1 = getCols(C1);
S2 = getCols(C2);
S3 = getCols(C3);

% Function handle to average repeated measurements at each load
avgByLoad = @(Load,y) localAvgByLoad(Load,y);

% Function handle for plotting all three cases with linear fits
plot3CasesWithFit = @(yname, ylab) localPlot3CasesWithFit(S1,S2,S3,avgByLoad,yname,ylab);

% Generate required plots (measurement vs applied load)
plot3CasesWithFit('F0' , 'F0 Reaction Force (lbf)');
plot3CasesWithFit('F1' , 'F1 Reaction Force (lbf)');
plot3CasesWithFit('F2' , 'F2 Reaction Force (lbf)');
plot3CasesWithFit('F3D', 'Internal Force F3D (lbf)');
plot3CasesWithFit('LVDT','LVDT Displacement (in)');

% Print regression parameters and R^2 values
localPrintFitStats(S1,S2,S3,avgByLoad,'F0');
localPrintFitStats(S1,S2,S3,avgByLoad,'F1');
localPrintFitStats(S1,S2,S3,avgByLoad,'F2');
localPrintFitStats(S1,S2,S3,avgByLoad,'F3D');
localPrintFitStats(S1,S2,S3,avgByLoad,'LVDT');

% Compute R^2 summary table for report
fields    = {'F0','F1','F2','F3D','LVDT'};
caseNames = {'Case 1','Case 2','Case 3'};
cases     = {S1,S2,S3};

R2_mat = zeros(numel(cases), numel(fields));

for i = 1:numel(cases)
    for j = 1:numel(fields)
        [L, ymean] = avgByLoad(cases{i}.Load, cases{i}.(fields{j}));
        p = polyfit(L, ymean, 1);
        R2_mat(i,j) = localR2(ymean, polyval(p,L));
    end
end

R2_Table = array2table(R2_mat, ...
    'VariableNames', fields, ...
    'RowNames', caseNames);

disp('----- R^2 Summary Table -----');
disp(R2_Table);

% ================= LOCAL FUNCTIONS =================

function [L, y_mean, y_std] = localAvgByLoad(Load, y)
    L = unique(Load);
    y_mean = zeros(size(L));
    y_std  = zeros(size(L));
    for i = 1:numel(L)
        idx = Load == L(i);
        y_mean(i) = mean(y(idx));
        y_std(i)  = std(y(idx));
    end
end

function localPlot3CasesWithFit(S1,S2,S3,avgByLoad,yfield,ylab)
    [L1, y1m, y1s] = avgByLoad(S1.Load, S1.(yfield));
    [L2, y2m, y2s] = avgByLoad(S2.Load, S2.(yfield));
    [L3, y3m, y3s] = avgByLoad(S3.Load, S3.(yfield));

    p1 = polyfit(L1, y1m, 1); y1fit = polyval(p1, L1); R21 = localR2(y1m, y1fit);
    p2 = polyfit(L2, y2m, 1); y2fit = polyval(p2, L2); R22 = localR2(y2m, y2fit);
    p3 = polyfit(L3, y3m, 1); y3fit = polyval(p3, L3); R23 = localR2(y3m, y3fit);

    figure; hold on; grid on;

    errorbar(L1, y1m, y1s, 'o', 'LineWidth',1.2);
    plot(L1, y1fit, '-', 'LineWidth',1.5);

    errorbar(L2, y2m, y2s, 's', 'LineWidth',1.2);
    plot(L2, y2fit, '-', 'LineWidth',1.5);

    errorbar(L3, y3m, y3s, '^', 'LineWidth',1.2);
    plot(L3, y3fit, '-', 'LineWidth',1.5);

    xlabel('Applied Load (lb)');
    ylabel(ylab);
    title(sprintf('%s vs Applied Load (3 Cases)', yfield));

    legend( ...
        sprintf('Case 1 data'), sprintf('Case 1 fit (R^2=%.3f)',R21), ...
        sprintf('Case 2 data'), sprintf('Case 2 fit (R^2=%.3f)',R22), ...
        sprintf('Case 3 data'), sprintf('Case 3 fit (R^2=%.3f)',R23), ...
        'Location','best');
end

function R2 = localR2(y, yfit)
    SSres = sum((y - yfit).^2);
    SStot = sum((y - mean(y)).^2);
    R2 = 1 - SSres/SStot;
end

function localPrintFitStats(S1,S2,S3,avgByLoad,yfield)
    [L1, y1m] = avgByLoad(S1.Load, S1.(yfield));
    [L2, y2m] = avgByLoad(S2.Load, S2.(yfield));
    [L3, y3m] = avgByLoad(S3.Load, S3.(yfield));

    p1 = polyfit(L1,y1m,1); R21 = localR2(y1m,polyval(p1,L1));
    p2 = polyfit(L2,y2m,1); R22 = localR2(y2m,polyval(p2,L2));
    p3 = polyfit(L3,y3m,1); R23 = localR2(y3m,polyval(p3,L3));

    fprintf('%s: Case1 slope=%.4g, intercept=%.4g, R^2=%.4f | ', yfield, p1(1), p1(2), R21);
    fprintf('Case2 slope=%.4g, intercept=%.4g, R^2=%.4f | ', yfield, p2(1), p2(2), R22);
    fprintf('Case3 slope=%.4g, intercept=%.4g, R^2=%.4f\n', yfield, p3(1), p3(2), R23);
end