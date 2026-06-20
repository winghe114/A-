%% Stochastic resonance k sensitivity experiment
% Sweep the input gain k with fixed a=b=1 and no extra artificial noise.

clear; clc; close all;

scriptDir = fileparts(mfilename('fullpath'));
q1Dir = fileparts(scriptDir);
rootDir = fileparts(q1Dir);
dataFile = fullfile(rootDir, 'data.xlsx');

fprintf('Running stochastic resonance k sensitivity experiment...\n');

[t, xRaw] = readSingleSourceData(dataFile);
t = t(:);
xRaw = xRaw(:);
y = xRaw - mean(xRaw, 'omitnan');
u = y / std(y, 0, 'omitnan');
fs = 1 / median(diff(t));
f0 = 2.0000016152054;
totalEnergy = sum(y .^ 2);

srParams = struct();
srParams.a = 1.0;
srParams.b = 1.0;
srParams.downsampleFactor = 2;
srParams.keepFraction = 0.10;
srParams.initialState = 0.05;
srParams.smoothingWindow = max(5, round(fs / srParams.downsampleFactor / f0 / 8));

kList = unique([linspace(0.005, 0.10, 28), linspace(0.12, 0.60, 49), linspace(0.65, 1.50, 18)])';
rows = cell(numel(kList), 11);
outputs = cell(numel(kList), 1);
times = cell(numel(kList), 1);

for i = 1:numel(kList)
    kGain = kList(i);
    [srTime, srOutput] = simulateSr(t, u, srParams, kGain);
    srProjectionFraction = coherentProjectionFraction(srTime, srOutput, f0);
    [srSnrDb, ~] = localFrequencySnr(srTime, srOutput, f0);

    rows(i, :) = {kGain, NaN, NaN, srProjectionFraction, srSnrDb, ...
        NaN, NaN, NaN, NaN, NaN, NaN};
    outputs{i} = srOutput;
    times{i} = srTime;
end

resultTable = cell2table(rows, 'VariableNames', {'k_gain', 'sr_frequency_Hz', ...
    'sr_coarse_frequency_Hz', 'sr_projection_fraction', 'sr_local_snr_dB', ...
    'original_projection_energy', 'original_residual_energy', ...
    'original_energy_fraction', 'difference_from_projection_Hz', ...
    'abs_difference_from_projection_Hz', 'combined_score'});

% Refine only a small set of promising gains. This keeps the sweep fast
% while still checking whether the best SR-looking gains give a better
% frequency when evaluated on the original signal.
srSnr = resultTable.sr_local_snr_dB;
srProj = resultTable.sr_projection_fraction;
zSnr = (srSnr - mean(srSnr, 'omitnan')) / std(srSnr, 0, 'omitnan');
zProj = (srProj - mean(srProj, 'omitnan')) / std(srProj, 0, 'omitnan');
quickScore = 0.5 * zSnr + 0.5 * zProj;
[~, sortIdx] = sort(quickScore, 'descend');
[~, idxMaxSnr] = max(resultTable.sr_local_snr_dB);
[~, idxMaxProj] = max(resultTable.sr_projection_fraction);
[~, idxNearestOld] = min(abs(resultTable.k_gain - 0.273333333333333));
candidateIdx = unique([sortIdx(1:min(16, numel(sortIdx))); idxMaxSnr; idxMaxProj; idxNearestOld], 'stable');

for c = 1:numel(candidateIdx)
    idx = candidateIdx(c);
    srTime = times{idx};
    srOutput = outputs{idx};
    [fSr, fSrCoarse, ~, ~] = projectionSearchLight(srTime, srOutput - mean(srOutput), f0, 0.35, 801);
    [~, fitOrig, resOrig] = fitSinusoidAtFrequency(t, y, fSr);
    origProjectionEnergy = sum(fitOrig .^ 2);
    origResidualEnergy = sum(resOrig .^ 2);
    origEnergyFraction = origProjectionEnergy / totalEnergy;
    score = 0.45 * resultTable.sr_projection_fraction(idx) + ...
        0.35 * resultTable.sr_local_snr_dB(idx) / 50 + ...
        0.20 * origEnergyFraction / 0.058294029412865;

    resultTable.sr_frequency_Hz(idx) = fSr;
    resultTable.sr_coarse_frequency_Hz(idx) = fSrCoarse;
    resultTable.original_projection_energy(idx) = origProjectionEnergy;
    resultTable.original_residual_energy(idx) = origResidualEnergy;
    resultTable.original_energy_fraction(idx) = origEnergyFraction;
    resultTable.difference_from_projection_Hz(idx) = fSr - f0;
    resultTable.abs_difference_from_projection_Hz(idx) = abs(fSr - f0);
    resultTable.combined_score(idx) = score;
end

writetable(resultTable, fullfile(scriptDir, 'sr_k_sensitivity_results.csv'));

makeKSensitivityFigure(scriptDir, resultTable, f0);
makeBestKPreviewFigure(scriptDir, t, u, resultTable, times, outputs, f0);
writeKSummary(scriptDir, resultTable);

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

function makeKSensitivityFigure(outDir, tbl, f0)
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 1180, 900]);
    tiledlayout(4, 1, 'TileSpacing', 'compact');
    nexttile;
    plot(tbl.k_gain, tbl.sr_local_snr_dB, 'Color', [0.75,0.12,0.18], 'LineWidth', 1.25);
    grid on;
    xlabel('Input gain k');
    ylabel('SR local SNR / dB');
    title('SR local SNR under different input gains');

    nexttile;
    plot(tbl.k_gain, tbl.sr_projection_fraction, 'Color', [0.1,0.45,0.72], 'LineWidth', 1.25);
    grid on;
    xlabel('Input gain k');
    ylabel('SR projection fraction');
    title('Projection fraction inside SR output');

    nexttile;
    validOrig = ~isnan(tbl.original_energy_fraction);
    plot(tbl.k_gain(validOrig), tbl.original_energy_fraction(validOrig), '-o', 'Color', [0.15,0.55,0.35], 'LineWidth', 1.25);
    grid on;
    xlabel('Input gain k');
    ylabel('Original projection fraction');
    title('SR-estimated frequency evaluated on original signal');

    nexttile;
    validDiff = ~isnan(tbl.difference_from_projection_Hz);
    plot(tbl.k_gain(validDiff), abs(tbl.difference_from_projection_Hz(validDiff)), '-o', 'Color', [0.42,0.23,0.62], 'LineWidth', 1.25);
    grid on;
    xlabel('Input gain k');
    ylabel('|Frequency difference| / Hz');
    title(sprintf('Difference from projection reference %.12f Hz', f0));
    exportgraphics(fig, fullfile(outDir, 'fig_sr_k_sensitivity.png'), 'Resolution', 220);
    close(fig);
end

function makeBestKPreviewFigure(outDir, t, u, tbl, times, outputs, f0)
    validScore = ~isnan(tbl.combined_score);
    [~, localScore] = max(tbl.combined_score(validScore));
    scoreIdxList = find(validScore);
    idxScore = scoreIdxList(localScore);
    [~, idxSnr] = max(tbl.sr_local_snr_dB);
    validDiff = ~isnan(tbl.abs_difference_from_projection_Hz);
    [~, localDiff] = min(tbl.abs_difference_from_projection_Hz(validDiff));
    diffIdxList = find(validDiff);
    idxDiff = diffIdxList(localDiff);
    indices = unique([idxScore, idxSnr, idxDiff], 'stable');

    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 1180, 360 * numel(indices)]);
    tiledlayout(numel(indices), 2, 'TileSpacing', 'compact');
    for n = 1:numel(indices)
        idx = indices(n);
        srTime = times{idx};
        srOutput = normalizeColumn(outputs{idx});
        tStart = srTime(1);
        tEnd = min(tStart + 8, srTime(end));
        showOriginal = t >= tStart & t <= tEnd;
        showSr = srTime >= tStart & srTime <= tEnd;

        nexttile;
        plot(t(showOriginal), u(showOriginal), 'Color', [0.58,0.58,0.58], 'LineWidth', 0.7);
        hold on;
        plot(srTime(showSr), srOutput(showSr), 'Color', [0.75,0.12,0.18], 'LineWidth', 1.1);
        grid on;
        xlabel('Time / s');
        ylabel('Normalized amplitude');
        title(sprintf('k=%.4f, f=%.9f Hz', tbl.k_gain(idx), tbl.sr_frequency_Hz(idx)));

        fs = 1 / median(diff(t));
        fsSr = 1 / median(diff(srTime));
        [fOrig, ampOrig] = amplitudeSpectrum(u, fs);
        [fSr, ampSr] = amplitudeSpectrum(srOutput, fsSr);
        ampOrig = normalizeSpectrumInBand(fOrig, ampOrig, 0, 8);
        ampSr = normalizeSpectrumInBand(fSr, ampSr, 0, 8);
        nexttile;
        plot(fOrig, ampOrig, 'Color', [0.58,0.58,0.58], 'LineWidth', 1.0);
        hold on;
        plot(fSr, ampSr, 'Color', [0.75,0.12,0.18], 'LineWidth', 1.2);
        xline(f0, '--', 'Color', [0.1,0.28,0.58], 'LineWidth', 1.1);
        xlim([0, 8]);
        grid on;
        xlabel('Frequency / Hz');
        ylabel('Normalized amplitude');
        title(sprintf('SNR %.2f dB, SR proj %.4f', tbl.sr_local_snr_dB(idx), tbl.sr_projection_fraction(idx)));
    end
    exportgraphics(fig, fullfile(outDir, 'fig_sr_k_best_cases.png'), 'Resolution', 220);
    close(fig);
end

function writeKSummary(outDir, tbl)
    validScore = ~isnan(tbl.combined_score);
    [~, localScore] = max(tbl.combined_score(validScore));
    scoreIdxList = find(validScore);
    idxScore = scoreIdxList(localScore);
    [~, idxSnr] = max(tbl.sr_local_snr_dB);
    [~, idxProj] = max(tbl.sr_projection_fraction);
    validOrig = ~isnan(tbl.original_energy_fraction);
    [~, localOrig] = max(tbl.original_energy_fraction(validOrig));
    origIdxList = find(validOrig);
    idxOrig = origIdxList(localOrig);
    validDiff = ~isnan(tbl.abs_difference_from_projection_Hz);
    [~, localDiff] = min(tbl.abs_difference_from_projection_Hz(validDiff));
    diffIdxList = find(validDiff);
    idxDiff = diffIdxList(localDiff);

    summary = table( ...
        ["combined_score"; "sr_local_snr"; "sr_projection_fraction"; "original_energy_fraction"; "min_frequency_difference"], ...
        [tbl.k_gain(idxScore); tbl.k_gain(idxSnr); tbl.k_gain(idxProj); tbl.k_gain(idxOrig); tbl.k_gain(idxDiff)], ...
        [tbl.sr_frequency_Hz(idxScore); tbl.sr_frequency_Hz(idxSnr); tbl.sr_frequency_Hz(idxProj); tbl.sr_frequency_Hz(idxOrig); tbl.sr_frequency_Hz(idxDiff)], ...
        [tbl.sr_projection_fraction(idxScore); tbl.sr_projection_fraction(idxSnr); tbl.sr_projection_fraction(idxProj); tbl.sr_projection_fraction(idxOrig); tbl.sr_projection_fraction(idxDiff)], ...
        [tbl.sr_local_snr_dB(idxScore); tbl.sr_local_snr_dB(idxSnr); tbl.sr_local_snr_dB(idxProj); tbl.sr_local_snr_dB(idxOrig); tbl.sr_local_snr_dB(idxDiff)], ...
        [tbl.original_energy_fraction(idxScore); tbl.original_energy_fraction(idxSnr); tbl.original_energy_fraction(idxProj); tbl.original_energy_fraction(idxOrig); tbl.original_energy_fraction(idxDiff)], ...
        [tbl.difference_from_projection_Hz(idxScore); tbl.difference_from_projection_Hz(idxSnr); tbl.difference_from_projection_Hz(idxProj); tbl.difference_from_projection_Hz(idxOrig); tbl.difference_from_projection_Hz(idxDiff)], ...
        'VariableNames', {'criterion', 'k_gain', 'sr_frequency_Hz', ...
        'sr_projection_fraction', 'sr_local_snr_dB', ...
        'original_energy_fraction', 'difference_from_projection_Hz'});
    writetable(summary, fullfile(outDir, 'sr_k_sensitivity_summary.csv'));
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
