clear; clc;

scriptDir = fileparts(mfilename('fullpath'));
projectDir = fileparts(scriptDir);
dataFile = fullfile(projectDir, 'data.xlsx');
outCsv = fullfile(scriptDir, 'frequency_domain_fft.csv');
outPng = fullfile(scriptDir, 'frequency_domain_fft.png');
outFig = fullfile(scriptDir, 'frequency_domain_fft.fig');

opts = detectImportOptions(dataFile, 'Sheet', 2, 'VariableNamingRule', 'preserve');
T = readtable(dataFile, opts);

t = T{:, 1};
x = T{:, 2};
valid = isfinite(t) & isfinite(x);
t = t(valid);
x = x(valid);

fs = 1 / median(diff(t));
x = x - mean(x);
N = numel(x);
Y = fft(x);
P2 = abs(Y / N);
P1 = P2(1:floor(N / 2) + 1);
P1(2:end-1) = 2 * P1(2:end-1);
f = fs * (0:floor(N / 2))' / N;

resultTable = table(f, P1, 'VariableNames', {'frequency_Hz', 'amplitude'});
writetable(resultTable, outCsv);

figure('Color', 'w', 'Position', [100, 100, 1100, 600]);
plot(f, P1, 'LineWidth', 1.1);
grid on;
title('x(t) 频域图（单边 FFT 幅值谱）');
xlabel('频率 / Hz');
ylabel('幅值');
xlim([0, fs / 2]);

exportgraphics(gcf, outPng, 'Resolution', 300);
savefig(gcf, outFig);
