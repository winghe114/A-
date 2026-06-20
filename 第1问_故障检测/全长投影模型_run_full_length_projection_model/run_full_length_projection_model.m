%% Full-length sinusoidal projection model for Question 1
% Generate figures and result tables for Section 5.1.2-5.1.3.

clear; clc; close all;

scriptDir = fileparts(mfilename('fullpath'));
q1Dir = fileparts(scriptDir);
rootDir = fileparts(q1Dir);
dataFile = fullfile(rootDir, 'data.xlsx');

fprintf('Running full-length projection model...\n');

[t, xRaw] = readSingleSourceData(dataFile);
t = t(:);
xRaw = xRaw(:);
y = xRaw - mean(xRaw, 'omitnan');
fs = 1 / median(diff(t));
n = numel(y);
duration = t(end) - t(1);
totalEnergy = sum(y .^ 2);

%% Full-length FFT coarse localization
[freq, amp, powerRaw] = fullLengthSpectrum(y, fs);
valid = freq > 0.05 & freq < min(49.5, fs / 2 - 0.1);
validIdx = find(valid);
[~, localMax] = max(powerRaw(validIdx));
fftIdx = validIdx(localMax);
fFft = freq(fftIdx);
[~, fitFft, resFft] = fitSinusoidAtFrequency(t, y, fFft);

%% Projection search near FFT peak
halfWidth = 0.04;
gridN = 6001;
fGrid = linspace(fFft - halfWidth, fFft + halfWidth, gridN)';
projEnergy = zeros(size(fGrid));
resEnergy = zeros(size(fGrid));
detectionRatio = zeros(size(fGrid));
for k = 1:numel(fGrid)
    [~, fit, res] = fitSinusoidAtFrequency(t, y, fGrid(k));
    projEnergy(k) = sum(fit .^ 2);
    resEnergy(k) = sum(res .^ 2);
    detectionRatio(k) = projEnergy(k) / max(resEnergy(k), eps);
end
[~, bestGridIdx] = max(projEnergy);
fCoarseProjection = fGrid(bestGridIdx);
gridStep = fGrid(2) - fGrid(1);
objective = @(f) -projectionEnergyAtFrequency(t, y, f);
opts = optimset('TolX', 1e-12, 'Display', 'off');
fProjection = fminbnd(objective, fCoarseProjection - 5*gridStep, fCoarseProjection + 5*gridStep, opts);
[thetaProjection, fitProjection, resProjection] = fitSinusoidAtFrequency(t, y, fProjection);
projectionEnergy = sum(fitProjection .^ 2);
projectionResidualEnergy = sum(resProjection .^ 2);
projectionDetectionRatio = projectionEnergy / projectionResidualEnergy;

%% Dense cos-sin correlation / GLRT verification
[fCorrelation, corrScore, fCorrelationCoarse, corrGrid, corrCurve] = ...
    denseCorrelationScan(t, y, fFft, halfWidth);
[~, fitCorr, resCorr] = fitSinusoidAtFrequency(t, y, fCorrelation);

%% Save CSV files
method = {'Full-length FFT coarse bin'; 'Full-length sinusoidal projection / GLRT'; 'Dense cos-sin correlation scan'};
frequency_Hz = [fFft; fProjection; fCorrelation];
coarse_frequency_Hz = [fFft; fCoarseProjection; fCorrelationCoarse];
projection_energy = [sum(fitFft .^ 2); projectionEnergy; sum(fitCorr .^ 2)];
residual_energy = [sum(resFft .^ 2); projectionResidualEnergy; sum(resCorr .^ 2)];
energy_fraction = projection_energy / totalEnergy;
detection_ratio = projection_energy ./ residual_energy;
difference_from_fft_Hz = frequency_Hz - fFft;
difference_from_projection_Hz = frequency_Hz - fProjection;
results = table(method, frequency_Hz, coarse_frequency_Hz, projection_energy, ...
    residual_energy, energy_fraction, detection_ratio, difference_from_fft_Hz, ...
    difference_from_projection_Hz);
writetable(results, fullfile(scriptDir, 'full_length_projection_results.csv'));

curveTable = table(fGrid, projEnergy, resEnergy, detectionRatio, ...
    'VariableNames', {'frequency_Hz', 'projection_energy_J', ...
    'residual_energy', 'detection_ratio'});
writetable(curveTable, fullfile(scriptDir, 'projection_search_curve.csv'));

corrTable = table(corrGrid, corrCurve, ...
    'VariableNames', {'frequency_Hz', 'correlation_GLRT_score'});
writetable(corrTable, fullfile(scriptDir, 'correlation_glrt_curve.csv'));

%% Figures
makeFftVsProjectionFigure(scriptDir, freq, amp, fFft, fProjection);
makeProjectionEnergyFigure(scriptDir, fGrid, projEnergy, fFft, fProjection);
makeCorrelationGlrtFigure(scriptDir, corrGrid, corrCurve, fFft, fProjection, fCorrelation);

fprintf('Finished. Outputs saved in %s\n', scriptDir);

%% Local functions
function [t, x] = readSingleSourceData(dataFile)
    try
        tbl = readtable(dataFile, 'Sheet', 1, 'VariableNamingRule', 'preserve');
        t = tbl{:, 1};
        x = tbl{:, 2};
    catch
        num = readmatrix(dataFile, 'Sheet', 1);
        t = num(:, 1);
        x = num(:, 2);
    end
    good = isfinite(t) & isfinite(x);
    t = t(good);
    x = x(good);
end

function [freq, amp, powerRaw] = fullLengthSpectrum(x, fs)
    n = numel(x);
    X = fft(x);
    freq = (0:floor(n / 2))' * fs / n;
    powerRaw = abs(X(1:numel(freq))).^2;
    amp = abs(X(1:numel(freq))) / n * 2;
    amp(1) = amp(1) / 2;
end

function [theta, fit, residual] = fitSinusoidAtFrequency(t, x, f)
    t = t(:);
    x = x(:) - mean(x, 'omitnan');
    H = [cos(2*pi*f*t), sin(2*pi*f*t)];
    theta = H \ x;
    fit = H * theta;
    residual = x - fit;
end

function e = projectionEnergyAtFrequency(t, x, f)
    [~, fit, ~] = fitSinusoidAtFrequency(t, x, f);
    e = sum(fit .^ 2);
end

function [fHat, bestScore, coarse, fGrid, scores] = denseCorrelationScan(t, x, center, halfWidth)
    fGrid = linspace(center - halfWidth, center + halfWidth, 4001)';
    scores = zeros(size(fGrid));
    x = x(:) - mean(x, 'omitnan');
    t = t(:);
    for i = 1:numel(fGrid)
        c = cos(2*pi*fGrid(i)*t);
        s = sin(2*pi*fGrid(i)*t);
        scores(i) = (x' * c)^2 / (c' * c) + (x' * s)^2 / (s' * s);
    end
    [bestScore, idx] = max(scores);
    coarse = fGrid(idx);
    objective = @(f) -projectionEnergyAtFrequency(t, x, f);
    opts = optimset('TolX', 1e-12, 'Display', 'off');
    fHat = fminbnd(objective, max(0.001, coarse - 0.005), coarse + 0.005, opts);
end

function makeFftVsProjectionFigure(outDir, freq, amp, fFft, fProjection)
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 1100, 680]);
    tiledlayout(2, 1, 'TileSpacing', 'compact');
    nexttile;
    plot(freq, amp, 'Color', [0.55,0.55,0.55], 'LineWidth', 0.9);
    hold on;
    xline(fFft, '--', 'Color', [0.1,0.45,0.72], 'LineWidth', 1.3);
    xline(fProjection, '-', 'Color', [0.75,0.12,0.18], 'LineWidth', 1.4);
    xlim([0, 8]);
    grid on;
    xlabel('Frequency / Hz');
    ylabel('Amplitude');
    title('Full-length FFT coarse localization and projection-refined frequency');
    legend({'Full-length FFT amplitude', 'FFT coarse bin', 'Projection estimate'}, 'Location', 'best');

    nexttile;
    plot(freq, amp, 'Color', [0.55,0.55,0.55], 'LineWidth', 0.9);
    hold on;
    xline(fFft, '--', 'Color', [0.1,0.45,0.72], 'LineWidth', 1.3);
    xline(fProjection, '-', 'Color', [0.75,0.12,0.18], 'LineWidth', 1.4);
    xlim([1.985, 2.015]);
    grid on;
    xlabel('Frequency / Hz');
    ylabel('Amplitude');
    title('Zoom near 2 Hz: FFT grid point versus continuous projection estimate');
    legend({'FFT amplitude', sprintf('FFT %.12f Hz', fFft), ...
        sprintf('Projection %.12f Hz', fProjection)}, 'Location', 'best');
    exportgraphics(fig, fullfile(outDir, 'fig_fft_vs_projection.png'), 'Resolution', 220);
    close(fig);
end

function makeProjectionEnergyFigure(outDir, fGrid, projEnergy, fFft, fProjection)
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 1100, 620]);
    plot(fGrid, projEnergy, 'Color', [0.75,0.12,0.18], 'LineWidth', 1.4);
    hold on;
    xline(fFft, '--', 'Color', [0.1,0.45,0.72], 'LineWidth', 1.2);
    xline(fProjection, '-', 'Color', [0.12,0.55,0.32], 'LineWidth', 1.3);
    grid on;
    xlabel('Frequency / Hz');
    ylabel('Projection energy J(f)');
    title('Full-length sinusoidal projection energy curve');
    legend({'J(f)', 'FFT coarse bin', 'Projection maximum'}, 'Location', 'best');
    exportgraphics(fig, fullfile(outDir, 'fig_projection_energy_curve.png'), 'Resolution', 220);
    close(fig);
end

function makeCorrelationGlrtFigure(outDir, fGrid, corrCurve, fFft, fProjection, fCorrelation)
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 1100, 620]);
    plot(fGrid, corrCurve, 'Color', [0.42,0.23,0.62], 'LineWidth', 1.4);
    hold on;
    xline(fFft, '--', 'Color', [0.1,0.45,0.72], 'LineWidth', 1.2);
    xline(fProjection, '-', 'Color', [0.75,0.12,0.18], 'LineWidth', 1.3);
    xline(fCorrelation, ':', 'Color', [0.12,0.55,0.32], 'LineWidth', 1.4);
    grid on;
    xlabel('Frequency / Hz');
    ylabel('Detection statistic');
    title('Correlation / GLRT statistic curve');
    legend({'Projection-to-residual energy ratio', 'FFT coarse bin', ...
        'Projection estimate', 'Correlation estimate'}, 'Location', 'best');
    exportgraphics(fig, fullfile(outDir, 'fig_correlation_glrt_curve.png'), 'Resolution', 220);
    close(fig);
end
