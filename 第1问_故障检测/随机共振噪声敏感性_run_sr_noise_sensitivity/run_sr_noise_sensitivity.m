%% Stochastic resonance noise sensitivity experiment
% This script tests whether adding controlled Gaussian noise improves the
% stochastic resonance response. It is an exploratory experiment and should
% not be used as the final estimator for Question 1.

clear; clc; close all;
rng(20260619);

scriptDir = fileparts(mfilename('fullpath'));
q1Dir = fileparts(scriptDir);
rootDir = fileparts(q1Dir);
dataFile = fullfile(rootDir, 'data.xlsx');

fprintf('Running stochastic resonance noise sensitivity experiment...\n');

[t, xRaw] = readSingleSourceData(dataFile);
t = t(:);
xRaw = xRaw(:);
y = xRaw - mean(xRaw, 'omitnan');
uBase = y / std(y, 0, 'omitnan');
dt = median(diff(t));
fs = 1 / dt;
f0 = 2.0000016152054;
totalEnergy = sum(y .^ 2);

srParams = struct();
srParams.a = 1.0;
srParams.b = 1.0;
srParams.downsampleFactor = 2;
srParams.keepFraction = 0.10;
srParams.initialState = 0.05;
srParams.smoothingWindow = max(5, round(fs / srParams.downsampleFactor / f0 / 8));

noiseSigmaList = [0, 0.05, 0.10, 0.20, 0.30, 0.50, 0.80, 1.00, 1.50, 2.00]';
gainList = linspace(0.02, 0.40, 25)';
rows = cell(numel(noiseSigmaList), 10);
bestOutputs = cell(numel(noiseSigmaList), 1);
bestTimes = cell(numel(noiseSigmaList), 1);

for i = 1:numel(noiseSigmaList)
    sigmaAdd = noiseSigmaList(i);
    addedNoise = sigmaAdd * randn(size(uBase));
    uNoisy = uBase + addedNoise;
    uNoisy = uNoisy - mean(uNoisy, 'omitnan');
    uNoisy = uNoisy / std(uNoisy, 0, 'omitnan');

    [bestGain, srOutput, srTime, srProjectionFraction, srSnrDb] = ...
        bestStochasticResonance(t, uNoisy, f0, srParams, gainList);
    [fSr, fSrCoarse, ~, ~] = projectionSearchLight(srTime, srOutput - mean(srOutput), f0, 0.35, 2001);
    [~, fitOrig, resOrig] = fitSinusoidAtFrequency(t, y, fSr);
    origProjectionEnergy = sum(fitOrig .^ 2);
    origResidualEnergy = sum(resOrig .^ 2);
    origEnergyFraction = origProjectionEnergy / totalEnergy;

    rows(i, :) = {sigmaAdd, bestGain, fSr, fSrCoarse, ...
        srProjectionFraction, srSnrDb, origProjectionEnergy, origResidualEnergy, ...
        origEnergyFraction, fSr - f0};
    bestOutputs{i} = srOutput;
    bestTimes{i} = srTime;
end

resultTable = cell2table(rows, 'VariableNames', {'added_noise_sigma', ...
    'best_gain', 'sr_frequency_Hz', 'sr_coarse_frequency_Hz', ...
    'sr_projection_fraction', 'sr_local_snr_dB', ...
    'original_projection_energy', 'original_residual_energy', ...
    'original_energy_fraction', 'difference_from_projection_Hz'});
writetable(resultTable, fullfile(scriptDir, 'sr_noise_sensitivity_results.csv'));

makeNoiseSensitivityFigure(scriptDir, resultTable, f0);
makeBestNoisePreviewFigure(scriptDir, t, uBase, resultTable, bestTimes, bestOutputs, f0);

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

function [bestGain, bestOutput, bestTime, bestProj, bestSnrDb] = bestStochasticResonance(t, xNorm, f0, p, gainList)
    bestScore = -Inf;
    bestGain = gainList(1);
    bestOutput = [];
    bestTime = [];
    bestProj = NaN;
    bestSnrDb = NaN;
    for i = 1:numel(gainList)
        [tout, y] = simulateSr(t, xNorm, p, gainList(i));
        proj = coherentProjectionFraction(tout, y, f0);
        [snrDb, ~] = localFrequencySnr(tout, y, f0);
        score = 0.55 * proj + 0.45 * snrDb / 50;
        if score > bestScore
            bestScore = score;
            bestGain = gainList(i);
            bestOutput = y;
            bestTime = tout;
            bestProj = proj;
            bestSnrDb = snrDb;
        end
    end
end

function [tout, yout] = simulateSr(t, xNorm, p, inputGain)
    ds = max(1, round(p.downsampleFactor));
    td = t(1:ds:end);
    u = xNorm(1:ds:end);
    if p.smoothingWindow > 1
        u = movmean(u, p.smoothingWindow);
    end
    dt = median(diff(td));
    y = p.initialState;
    keepStart = max(1, floor(numel(td) * p.keepFraction));
    yout = zeros(numel(td)-keepStart, 1);
    tout = zeros(numel(td)-keepStart, 1);
    outIdx = 0;
    for k = 1:numel(td)-1
        inputValue = inputGain * u(k);
        k1 = srDeriv(y, inputValue, p);
        k2 = srDeriv(y + 0.5*dt*k1, inputValue, p);
        k3 = srDeriv(y + 0.5*dt*k2, inputValue, p);
        k4 = srDeriv(y + dt*k3, inputValue, p);
        y = y + dt * (k1 + 2*k2 + 2*k3 + k4) / 6;
        if abs(y) > 10
            y = sign(y) * 10;
        end
        if k >= keepStart
            outIdx = outIdx + 1;
            tout(outIdx) = td(k);
            yout(outIdx) = y;
        end
    end
    tout = tout(1:outIdx);
    yout = yout(1:outIdx);
    yout = yout - mean(yout, 'omitnan');
end

function dy = srDeriv(y, inputValue, p)
    dy = p.a * y - p.b * y^3 + inputValue;
end

function frac = coherentProjectionFraction(t, y, f)
    y = y(:) - mean(y);
    t = t(:);
    H = [cos(2*pi*f*t), sin(2*pi*f*t)];
    theta = H \ y;
    fit = H * theta;
    frac = sum(fit .^ 2) / max(sum(y .^ 2), eps);
end

function [snrDb, peakAmp] = localFrequencySnr(t, y, f0)
    y = y(:) - mean(y);
    fs = 1 / median(diff(t));
    [freq, amp] = amplitudeSpectrum(y, fs);
    [~, idx] = min(abs(freq - f0));
    lo = max(2, idx - 30);
    hi = min(numel(amp), idx + 30);
    localAmp = amp(lo:hi);
    mask = true(size(localAmp));
    center = idx - lo + 1;
    mask(max(1, center-2):min(numel(mask), center+2)) = false;
    noiseFloor = median(localAmp(mask));
    peakAmp = amp(idx);
    snrDb = 20 * log10(peakAmp / max(noiseFloor, eps));
end

function [theta, fit, residual] = fitSinusoidAtFrequency(t, x, f)
    t = t(:);
    x = x(:) - mean(x);
    H = [cos(2*pi*f*t), sin(2*pi*f*t)];
    theta = H \ x;
    fit = H * theta;
    residual = x - fit;
end

function [fHat, coarseFreq, projEnergy, residualEnergy] = projectionSearchLight(t, x, center, halfWidth, gridN)
    t = t(:);
    x = x(:) - mean(x);
    dt = median(diff(t));
    fs = 1 / dt;
    lo = max(0.001, center - halfWidth);
    hi = min(fs / 2 - 0.001, center + halfWidth);
    fGrid = linspace(lo, hi, gridN)';
    eGrid = zeros(size(fGrid));
    for k = 1:numel(fGrid)
        [~, fit, ~] = fitSinusoidAtFrequency(t, x, fGrid(k));
        eGrid(k) = sum(fit .^ 2);
    end
    [~, bestIdx] = max(eGrid);
    coarseFreq = fGrid(bestIdx);
    gridStep = fGrid(2) - fGrid(1);
    refineLo = max(lo, coarseFreq - 5 * gridStep);
    refineHi = min(hi, coarseFreq + 5 * gridStep);
    objective = @(f) -projectionEnergy(t, x, f);
    opts = optimset('TolX', 1e-12, 'Display', 'off');
    fHat = fminbnd(objective, refineLo, refineHi, opts);
    projEnergy = -objective(fHat);
    [~, ~, residual] = fitSinusoidAtFrequency(t, x, fHat);
    residualEnergy = sum(residual .^ 2);
end

function e = projectionEnergy(t, x, f)
    [~, fit, ~] = fitSinusoidAtFrequency(t, x, f);
    e = sum(fit .^ 2);
end

function [freq, amp] = amplitudeSpectrum(x, fs)
    x = x(:) - mean(x(:), 'omitnan');
    n = numel(x);
    win = hannLocal(n);
    X = fft(x .* win);
    freq = (0:floor(n/2))' * fs / n;
    amp = abs(X(1:numel(freq))) / sum(win) * 2;
end

function w = hannLocal(n)
    k = (0:n-1)';
    w = 0.5 - 0.5 * cos(2*pi*k/max(n-1, 1));
end

function makeNoiseSensitivityFigure(outDir, tbl, f0)
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 1180, 820]);
    tiledlayout(3, 1, 'TileSpacing', 'compact');
    nexttile;
    plot(tbl.added_noise_sigma, tbl.sr_local_snr_dB, '-o', 'Color', [0.75,0.12,0.18], 'LineWidth', 1.4);
    grid on;
    xlabel('Added Gaussian noise standard deviation');
    ylabel('SR local SNR / dB');
    title('Stochastic resonance local SNR under added noise');

    nexttile;
    plot(tbl.added_noise_sigma, tbl.original_energy_fraction, '-o', 'Color', [0.1,0.45,0.72], 'LineWidth', 1.4);
    grid on;
    xlabel('Added Gaussian noise standard deviation');
    ylabel('Projection fraction on original signal');
    title('Quality of SR-estimated frequency evaluated on original signal');

    nexttile;
    plot(tbl.added_noise_sigma, abs(tbl.difference_from_projection_Hz), '-o', 'Color', [0.15,0.55,0.35], 'LineWidth', 1.4);
    yline(0, '--', 'Color', [0.3,0.3,0.3]);
    grid on;
    xlabel('Added Gaussian noise standard deviation');
    ylabel('|Frequency difference| / Hz');
    title(sprintf('Frequency difference from projection reference %.12f Hz', f0));
    exportgraphics(fig, fullfile(outDir, 'fig_sr_noise_sensitivity.png'), 'Resolution', 220);
    close(fig);
end

function makeBestNoisePreviewFigure(outDir, t, uBase, tbl, bestTimes, bestOutputs, f0)
    [~, bestIdx] = max(tbl.sr_local_snr_dB);
    srTime = bestTimes{bestIdx};
    srOutput = bestOutputs{bestIdx};
    srOutput = normalizeColumn(srOutput);
    tStart = srTime(1);
    tEnd = min(tStart + 12, srTime(end));
    showOriginal = t >= tStart & t <= tEnd;
    showSr = srTime >= tStart & srTime <= tEnd;

    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 1180, 680]);
    tiledlayout(2, 1, 'TileSpacing', 'compact');
    nexttile;
    plot(t(showOriginal), uBase(showOriginal), 'Color', [0.58,0.58,0.58], 'LineWidth', 0.75);
    hold on;
    plot(srTime(showSr), srOutput(showSr), 'Color', [0.75,0.12,0.18], 'LineWidth', 1.15);
    grid on;
    xlabel('Time / s');
    ylabel('Normalized amplitude');
    title(sprintf('Best added-noise case, sigma = %.2f', tbl.added_noise_sigma(bestIdx)));
    legend({'Original standardized signal', 'SR output'}, 'Location', 'best');

    fs = 1 / median(diff(t));
    fsSr = 1 / median(diff(srTime));
    [fOrig, ampOrig] = amplitudeSpectrum(uBase, fs);
    [fSr, ampSr] = amplitudeSpectrum(srOutput, fsSr);
    ampOrig = normalizeSpectrumInBand(fOrig, ampOrig, 0, 8);
    ampSr = normalizeSpectrumInBand(fSr, ampSr, 0, 8);
    nexttile;
    plot(fOrig, ampOrig, 'Color', [0.58,0.58,0.58], 'LineWidth', 1.0);
    hold on;
    plot(fSr, ampSr, 'Color', [0.75,0.12,0.18], 'LineWidth', 1.2);
    xline(f0, '--', 'Color', [0.1,0.28,0.58], 'LineWidth', 1.2);
    xlim([0, 8]);
    grid on;
    xlabel('Frequency / Hz');
    ylabel('Normalized amplitude');
    title('Spectrum of the best added-noise SR case');
    legend({'Original spectrum', 'SR output spectrum', 'Projection reference'}, 'Location', 'best');
    exportgraphics(fig, fullfile(outDir, 'fig_sr_best_added_noise_preview.png'), 'Resolution', 220);
    close(fig);
end

function y = normalizeColumn(x)
    y = x(:) - mean(x(:), 'omitnan');
    scale = std(y, 0, 'omitnan');
    if ~isfinite(scale) || scale <= eps
        scale = 1;
    end
    y = y / scale;
end

function ampNorm = normalizeSpectrumInBand(freq, amp, lo, hi)
    ampNorm = amp;
    band = freq >= lo & freq <= hi;
    scale = max(ampNorm(band));
    if ~isfinite(scale) || scale <= eps
        scale = max(ampNorm);
    end
    if ~isfinite(scale) || scale <= eps
        scale = 1;
    end
    ampNorm = ampNorm / scale;
end
