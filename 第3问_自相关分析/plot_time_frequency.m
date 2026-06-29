clear; clc;

scriptDir = fileparts(mfilename('fullpath'));
projectDir = fileparts(scriptDir);
dataFile = fullfile(projectDir, 'data.xlsx');
outPng = fullfile(scriptDir, 'time_frequency_spectrogram.png');
outFig = fullfile(scriptDir, 'time_frequency_spectrogram.fig');

opts = detectImportOptions(dataFile, 'Sheet', 2, 'VariableNamingRule', 'preserve');
T = readtable(dataFile, opts);

t = T{:, 1};
x = T{:, 2};
valid = isfinite(t) & isfinite(x);
t = t(valid);
x = x(valid);

fs = 1 / median(diff(t));
windowLength = 1024;
overlapLength = round(0.75 * windowLength);
nfft = 2048;

[s, f, tt] = spectrogram(x, hamming(windowLength), overlapLength, nfft, fs);
powerDb = 20 * log10(abs(s) + eps);
timeAxis = tt + t(1);

figure('Color', 'w', 'Position', [100, 100, 1000, 600]);
imagesc(timeAxis, f, powerDb);
axis xy;
title('x(t) 时间-频率图（横轴为时间）');
xlabel('时间 / s');
ylabel('频率 / Hz');
colormap turbo;
h = colorbar;
ylabel(h, '幅值 / dB');
xlim([min(t), max(t)]);
ylim([0, fs / 2]);

exportgraphics(gcf, outPng, 'Resolution', 300);
savefig(gcf, outFig);
