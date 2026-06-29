clear; clc;

scriptDir = fileparts(mfilename('fullpath'));
projectDir = fileparts(scriptDir);
dataFile = fullfile(projectDir, 'data.xlsx');
outCsv = fullfile(scriptDir, 'sliding_autocorrelation_lag12.csv');
outPng = fullfile(scriptDir, 'sliding_autocorrelation_lag12.png');
outFig = fullfile(scriptDir, 'sliding_autocorrelation_lag12.fig');

opts = detectImportOptions(dataFile, 'Sheet', 2, 'VariableNamingRule', 'preserve');
T = readtable(dataFile, opts);

t = T{:, 1};
x = T{:, 2};
valid = isfinite(t) & isfinite(x);
t = t(valid);
x = x(valid);

windowLength = 501;
halfWindow = floor(windowLength / 2);
centerIdx = (1 + halfWindow):(numel(x) - halfWindow);
time = t(centerIdx);
lag1 = zeros(numel(centerIdx), 1);
lag2 = zeros(numel(centerIdx), 1);

for i = 1:numel(centerIdx)
    segment = x(centerIdx(i) - halfWindow:centerIdx(i) + halfWindow);
    segment = segment - mean(segment);
    denominator = sum(segment .^ 2);
    lag1(i) = sum(segment(1:end-1) .* segment(2:end)) / denominator;
    lag2(i) = sum(segment(1:end-2) .* segment(3:end)) / denominator;
end

resultTable = table(time, lag1, lag2, 'VariableNames', {'time', 'lag1_autocorrelation', 'lag2_autocorrelation'});
writetable(resultTable, outCsv);

figure('Color', 'w', 'Position', [100, 100, 1100, 600]);
plot(time, lag1, 'LineWidth', 1.2); hold on;
plot(time, lag2, 'LineWidth', 1.2);
grid on;
title('x(t) 滑动窗口一阶/二阶自相关');
xlabel('时间 / s');
ylabel('自相关值');
legend('一阶自相关', '二阶自相关', 'Location', 'best');
xlim([min(time), max(time)]);

exportgraphics(gcf, outPng, 'Resolution', 300);
savefig(gcf, outFig);
