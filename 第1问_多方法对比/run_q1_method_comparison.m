%% Question 1 multi-method comparison
% This script compares several methods that may outperform plain FFT for
% weak single-frequency fault detection.

clear; clc; close all;

scriptDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(scriptDir);
dataFile = fullfile(rootDir, 'data.xlsx');

fprintf('Question 1 multi-method comparison started...\n');

[t, xRaw] = readSingleSourceData(dataFile);
t = t(:);
xRaw = xRaw(:);
x = xRaw - mean(xRaw, 'omitnan');
dt = median(diff(t));
fs = 1 / dt;
duration = t(end) - t(1);
dfTheory = 1 / duration;
n = numel(x);
totalEnergy = sum(x .^ 2);

%% Plain full-length FFT baseline
[freq, amp, powerRaw] = fullLengthSpectrum(x, fs);
valid = freq > 0.05 & freq < min(49.5, fs / 2 - 0.1);
[~, validMax] = max(powerRaw(valid));
validIdx = find(valid);
fftIdx = validIdx(validMax);
fFftBin = freq(fftIdx);
[~, fitFft, resFft] = fitSinusoidAtFrequency(t, x, fFftBin);
fftEnergy = sum(fitFft .^ 2);

rows = {};
rows(end+1, :) = makeRow('Plain full-length FFT bin', fFftBin, fFftBin, ...
    fftEnergy, sum(resFft .^ 2), totalEnergy, 'FFT bin only; has grid bias.');

%% Quadratic interpolated FFT
fQuad = quadraticInterpolatedFftFrequency(freq, powerRaw, fftIdx);
[~, fitQuad, resQuad] = fitSinusoidAtFrequency(t, x, fQuad);
rows(end+1, :) = makeRow('Quadratic interpolated FFT', fQuad, fFftBin, ...
    sum(fitQuad .^ 2), sum(resQuad .^ 2), totalEnergy, ...
    'Low-cost correction of FFT grid bias.');

%% Zero-padded FFT
padFactor = 16;
[fZp, fZpCoarse] = zeroPaddedFftFrequency(x, fs, padFactor);
[~, fitZp, resZp] = fitSinusoidAtFrequency(t, x, fZp);
rows(end+1, :) = makeRow(sprintf('Zero-padded FFT x%d', padFactor), fZp, fZpCoarse, ...
    sum(fitZp .^ 2), sum(resZp .^ 2), totalEnergy, ...
    'Densifies the spectrum display but does not add information.');

%% Full-length projection / GLRT
[fProj, fProjCoarse, eProj, rProj] = projectionSearch(t, x, fFftBin, 0.04);
rows(end+1, :) = makeRow('Full-length projection / GLRT', fProj, fProjCoarse, ...
    eProj, rProj, totalEnergy, ...
    'Best main method; continuous-frequency coherent detection.');

%% Dense correlation scan
[fCorr, corrScore, corrCoarse] = denseCorrelationScan(t, x, fFftBin, 0.04);
[~, fitCorr, resCorr] = fitSinusoidAtFrequency(t, x, fCorr);
rows(end+1, :) = makeRow('Dense cos-sin correlation scan', fCorr, corrCoarse, ...
    sum(fitCorr .^ 2), sum(resCorr .^ 2), totalEnergy, ...
    sprintf('Equivalent to matched filtering; peak score %.6g.', corrScore));

%% AR prewhitening + projection
arOrder = 24;
[xWhite, arCoeff] = arPrewhiten(x, arOrder);
[fWhite, fWhiteCoarse, eWhite, rWhite] = projectionSearch(t(arOrder+1:end), xWhite, fFftBin, 0.04);
[~, fitWhiteOrig, resWhiteOrig] = fitSinusoidAtFrequency(t, x, fWhite);
rows(end+1, :) = makeRow(sprintf('AR(%d) prewhitened projection', arOrder), fWhite, fWhiteCoarse, ...
    sum(fitWhiteOrig .^ 2), sum(resWhiteOrig .^ 2), totalEnergy, ...
    'Tests robustness after reducing colored noise correlation.');

%% SSA denoising + projection
ssaWindow = 400;
ssaRank = 2;
xSsa = ssaDenoise(x, ssaWindow, ssaRank);
[fSsa, fSsaCoarse, eSsa, rSsa] = projectionSearch(t, xSsa - mean(xSsa), fFftBin, 0.04);
[~, fitSsaOrig, resSsaOrig] = fitSinusoidAtFrequency(t, x, fSsa);
rows(end+1, :) = makeRow(sprintf('SSA denoise rank-%d projection', ssaRank), fSsa, fSsaCoarse, ...
    sum(fitSsaOrig .^ 2), sum(resSsaOrig .^ 2), totalEnergy, ...
    sprintf('Data-adaptive low-rank denoising, window %d.', ssaWindow));

%% Stochastic resonance + projection
srParams = struct();
srParams.a = 1.0;
srParams.b = 1.0;
srParams.downsampleFactor = 2;
srParams.keepFraction = 0.10;
srParams.initialState = 0.05;
srParams.smoothingWindow = max(5, round(fs / srParams.downsampleFactor / fProj / 8));
gainList = linspace(0.02, 0.40, 25)';
[bestGain, srOutput, srTime, srProjectionFraction, srSnrDb] = ...
    bestStochasticResonance(t, x / std(x), fProj, srParams, gainList);
[fSr, fSrCoarse, eSr, rSr] = projectionSearch(srTime, srOutput - mean(srOutput), fProj, 0.35);
[~, fitSrOrig, resSrOrig] = fitSinusoidAtFrequency(t, x, fSr);
rows(end+1, :) = makeRow('Adaptive stochastic resonance + projection', fSr, fSrCoarse, ...
    sum(fitSrOrig .^ 2), sum(resSrOrig .^ 2), totalEnergy, ...
    sprintf('Nonlinear enhancement; best gain %.4g, output SNR %.3f dB.', bestGain, srSnrDb));

%% Multi-stage autocorrelation enhancement + projection
% Autocorrelation preserves the period of a sinusoid while suppressing
% uncorrelated random components in expectation. Repeating the operation can
% make the periodic structure more visible, but it also shortens the usable
% lag record and may introduce envelope bias, so the final comparison is
% still evaluated on the original time series.
acfMaxOrder = 3;
acfMaxLagSeconds = min(200, duration / 2);
[acfTimes, acfSignals] = multiAutocorrelationEnhancement(x, fs, acfMaxOrder, acfMaxLagSeconds);
acfMetricRows = {};
for acfOrder = 1:acfMaxOrder
    tau = acfTimes{acfOrder};
    yAcf = acfSignals{acfOrder};
    [fAcf, fAcfCoarse, eAcf, ~] = projectionSearchLight(tau, yAcf, fFftBin, 0.04, 2001);
    acfSelfFraction = eAcf / max(sum((yAcf - mean(yAcf)) .^ 2), eps);
    [acfSnrDb, ~] = localFrequencySnr(tau, yAcf, fAcf);
    [~, fitAcfOrig, resAcfOrig] = fitSinusoidAtFrequency(t, x, fAcf);
    rows(end+1, :) = makeRow(sprintf('Autocorrelation order-%d + projection', acfOrder), ...
        fAcf, fAcfCoarse, sum(fitAcfOrig .^ 2), sum(resAcfOrig .^ 2), totalEnergy, ...
        sprintf('Time-domain periodicity enhancement; autocorr-domain fraction %.6g, local SNR %.3f dB.', ...
        acfSelfFraction, acfSnrDb));
    acfMetricRows(end+1, :) = {acfOrder, fAcf, fAcfCoarse, acfSelfFraction, ...
        acfSnrDb, numel(yAcf), max(tau)};
end

%% Optional MUSIC if Signal Processing Toolbox is available
if exist('pmusic', 'file') == 2
    try
        [fMusic, fMusicCoarse] = musicEstimate(x, fs, fFftBin, 0.20);
        [~, fitMusic, resMusic] = fitSinusoidAtFrequency(t, x, fMusic);
        rows(end+1, :) = makeRow('MUSIC pseudospectrum', fMusic, fMusicCoarse, ...
            sum(fitMusic .^ 2), sum(resMusic .^ 2), totalEnergy, ...
            'Subspace high-resolution spectrum; toolbox dependent.');
    catch err
        rows(end+1, :) = makeRow('MUSIC pseudospectrum', NaN, NaN, NaN, NaN, totalEnergy, ...
            ['Skipped after runtime error: ', err.message]);
    end
else
    rows(end+1, :) = makeRow('MUSIC pseudospectrum', NaN, NaN, NaN, NaN, totalEnergy, ...
        'Skipped: pmusic is not available in this MATLAB environment.');
end

%% Save results
resultTable = cell2table(rows, 'VariableNames', {'method', 'frequency_Hz', ...
    'coarse_frequency_Hz', 'projection_energy', 'residual_energy', ...
    'energy_fraction', 'detection_ratio', 'difference_from_FFT_bin_Hz', ...
    'difference_from_projection_Hz', 'comment'});
fReference = resultTable.frequency_Hz(strcmp(resultTable.method, 'Full-length projection / GLRT'));
resultTable.difference_from_FFT_bin_Hz = resultTable.frequency_Hz - fFftBin;
resultTable.difference_from_projection_Hz = resultTable.frequency_Hz - fReference;
writetable(resultTable, fullfile(scriptDir, 'q1_method_comparison_results.csv'));

acfMetricTable = cell2table(acfMetricRows, 'VariableNames', {'order', ...
    'frequency_Hz', 'coarse_frequency_Hz', 'autocorrelation_projection_fraction', ...
    'local_snr_dB', 'sample_count', 'max_lag_s'});
writetable(acfMetricTable, fullfile(scriptDir, 'q1_autocorrelation_enhancement_metrics.csv'));

makeComparisonFigures(scriptDir, t, x, freq, amp, fFftBin, fProj, resultTable, ...
    xSsa, srTime, srOutput, acfTimes, acfSignals);
writeAnalysisReport(scriptDir, resultTable, fs, duration, dfTheory, fFftBin, fProj, ...
    bestGain, srProjectionFraction, srSnrDb);

fprintf('Comparison finished. Results saved in %s\n', scriptDir);

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

function row = makeRow(method, f, coarse, energy, residualEnergy, totalEnergy, comment)
    if isnan(f)
        energyFraction = NaN;
        detectionRatio = NaN;
        diffFft = NaN;
    else
        energyFraction = energy / totalEnergy;
        detectionRatio = energy / max(totalEnergy - energy, eps);
        diffFft = NaN;
    end
    row = {method, f, coarse, energy, residualEnergy, energyFraction, detectionRatio, diffFft, NaN, comment};
end

function [freq, amp, powerRaw] = fullLengthSpectrum(x, fs)
    n = numel(x);
    X = fft(x);
    freq = (0:floor(n / 2))' * fs / n;
    powerRaw = abs(X(1:numel(freq))).^2;
    amp = abs(X(1:numel(freq))) / n * 2;
    amp(1) = amp(1) / 2;
end

function f = quadraticInterpolatedFftFrequency(freq, powerRaw, idx)
    if idx <= 1 || idx >= numel(powerRaw)
        f = freq(idx);
        return;
    end
    y1 = log(max(powerRaw(idx-1), eps));
    y2 = log(max(powerRaw(idx), eps));
    y3 = log(max(powerRaw(idx+1), eps));
    denom = y1 - 2*y2 + y3;
    if abs(denom) < eps
        delta = 0;
    else
        delta = 0.5 * (y1 - y3) / denom;
    end
    df = freq(2) - freq(1);
    f = freq(idx) + delta * df;
end

function [fZp, coarse] = zeroPaddedFftFrequency(x, fs, padFactor)
    n = numel(x);
    nfft = 2^nextpow2(n * padFactor);
    X = fft(x, nfft);
    freq = (0:floor(nfft / 2))' * fs / nfft;
    p = abs(X(1:numel(freq))).^2;
    valid = freq > 0.05 & freq < min(49.5, fs / 2 - 0.1);
    validIdx = find(valid);
    [~, local] = max(p(validIdx));
    idx = validIdx(local);
    coarse = freq(idx);
    fZp = quadraticInterpolatedFftFrequency(freq, p, idx);
end

function [theta, fit, residual] = fitSinusoidAtFrequency(t, x, f)
    t = t(:);
    x = x(:) - mean(x);
    H = [cos(2*pi*f*t), sin(2*pi*f*t)];
    theta = H \ x;
    fit = H * theta;
    residual = x - fit;
end

function [fHat, coarseFreq, projEnergy, residualEnergy] = projectionSearch(t, x, center, halfWidth)
    t = t(:);
    x = x(:) - mean(x);
    dt = median(diff(t));
    fs = 1 / dt;
    lo = max(0.001, center - halfWidth);
    hi = min(fs / 2 - 0.001, center + halfWidth);
    gridN = 6001;
    fGrid = linspace(lo, hi, gridN)';
    eGrid = zeros(size(fGrid));
    for k = 1:numel(fGrid)
        eGrid(k) = projectionEnergy(t, x, fGrid(k));
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
        eGrid(k) = projectionEnergy(t, x, fGrid(k));
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

function [fHat, bestScore, coarse] = denseCorrelationScan(t, x, center, halfWidth)
    fGrid = linspace(center - halfWidth, center + halfWidth, 4001)';
    scores = zeros(size(fGrid));
    x = x(:) - mean(x);
    t = t(:);
    for i = 1:numel(fGrid)
        c = cos(2*pi*fGrid(i)*t);
        s = sin(2*pi*fGrid(i)*t);
        scores(i) = (x' * c)^2 / (c' * c) + (x' * s)^2 / (s' * s);
    end
    [bestScore, idx] = max(scores);
    coarse = fGrid(idx);
    localObjective = @(f) -projectionEnergy(t, x, f);
    opts = optimset('TolX', 1e-12, 'Display', 'off');
    fHat = fminbnd(localObjective, max(0.001, coarse - 0.005), coarse + 0.005, opts);
end

function [xWhite, arCoeff] = arPrewhiten(x, order)
    x = x(:) - mean(x);
    r = zeros(order+1, 1);
    for k = 0:order
        r(k+1) = sum(x(1:end-k) .* x(1+k:end)) / numel(x);
    end
    R = toeplitz(r(1:order));
    rhs = r(2:order+1);
    arCoeff = R \ rhs;
    pred = zeros(numel(x)-order, 1);
    for i = order+1:numel(x)
        past = x(i-1:-1:i-order);
        pred(i-order) = arCoeff' * past;
    end
    xWhite = x(order+1:end) - pred;
    xWhite = xWhite - mean(xWhite);
end

function xRecon = ssaDenoise(x, window, rankKeep)
    x = x(:) - mean(x);
    n = numel(x);
    L = min(window, floor(n / 2));
    K = n - L + 1;
    X = zeros(L, K);
    for k = 1:K
        X(:, k) = x(k:k+L-1);
    end
    [U, S, V] = svd(X, 'econ');
    rankKeep = min(rankKeep, size(S, 1));
    Xr = U(:, 1:rankKeep) * S(1:rankKeep, 1:rankKeep) * V(:, 1:rankKeep)';
    xRecon = diagonalAverage(Xr);
    xRecon = xRecon(:) - mean(xRecon);
end

function y = diagonalAverage(X)
    [L, K] = size(X);
    n = L + K - 1;
    y = zeros(n, 1);
    counts = zeros(n, 1);
    for i = 1:L
        for j = 1:K
            y(i+j-1) = y(i+j-1) + X(i, j);
            counts(i+j-1) = counts(i+j-1) + 1;
        end
    end
    y = y ./ counts;
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
    [freq, amp] = windowedAmplitude(y, fs);
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

function [freq, amp] = windowedAmplitude(x, fs)
    n = numel(x);
    win = hannLocal(n);
    X = fft(x(:) .* win);
    freq = (0:floor(n/2))' * fs / n;
    amp = abs(X(1:numel(freq))) / sum(win) * 2;
end

function w = hannLocal(n)
    k = (0:n-1)';
    w = 0.5 - 0.5 * cos(2*pi*k/max(n-1, 1));
end

function [fMusic, coarse] = musicEstimate(x, fs, center, halfWidth)
    x = x(:) - mean(x);
    [pxx, f] = pmusic(x, 2, 8192, fs);
    valid = f > center - halfWidth & f < center + halfWidth;
    idxValid = find(valid);
    [~, local] = max(pxx(idxValid));
    coarse = f(idxValid(local));
    fMusic = coarse;
end

function [acfTimes, acfSignals] = multiAutocorrelationEnhancement(x, fs, maxOrder, maxLagSeconds)
    acfTimes = cell(maxOrder, 1);
    acfSignals = cell(maxOrder, 1);
    current = normalizeColumn(x);
    for order = 1:maxOrder
        r = positiveAutocorrelationFft(current);
        if numel(r) > 1
            r = r(2:end);  % remove the dominant zero-lag spike
        end
        keepCount = min(numel(r), max(10, floor(maxLagSeconds * fs)));
        r = r(1:keepCount);
        tau = (1:keepCount)' / fs;
        r = detrend(r, 'linear');
        r = normalizeColumn(r);
        acfTimes{order} = tau;
        acfSignals{order} = r;
        current = r;
    end
end

function r = positiveAutocorrelationFft(x)
    x = x(:) - mean(x);
    n = numel(x);
    nfft = 2 ^ nextpow2(2 * n - 1);
    X = fft(x, nfft);
    r = real(ifft(X .* conj(X)));
    r = r(1:n) / n;  % biased estimate is more stable at large lags
end

function y = normalizeColumn(x)
    y = x(:) - mean(x(:), 'omitnan');
    scale = std(y, 'omitnan');
    if ~isfinite(scale) || scale <= eps
        scale = 1;
    end
    y = y / scale;
end
function makeComparisonFigures(outDir, t, x, freq, amp, fFftBin, fProj, resultTable, xSsa, srTime, srOutput, acfTimes, acfSignals)
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 1100, 620]);
    plot(freq, amp, 'Color', [0.1, 0.28, 0.58], 'LineWidth', 1.0);
    hold on;
    xline(fFftBin, 'Color', [0.5,0.5,0.5], 'LineStyle', '--', 'LineWidth', 1.0);
    xline(fProj, 'r--', 'LineWidth', 1.4);
    xlim([0, 10]);
    grid on;
    xlabel('Frequency / Hz');
    ylabel('Amplitude');
    title('Full-length FFT baseline and projection-refined frequency');
    legend({'FFT amplitude', 'FFT bin', 'Projection frequency'}, 'Location', 'best');
    exportgraphics(fig, fullfile(outDir, 'fig_baseline_fft_projection.png'), 'Resolution', 180);
    close(fig);

    validRows = ~isnan(resultTable.frequency_Hz);
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 1100, 620]);
    bar(categorical(resultTable.method(validRows)), resultTable.frequency_Hz(validRows));
    yline(fProj, 'r--', 'LineWidth', 1.2);
    grid on;
    ylabel('Estimated frequency / Hz');
    title('Frequency estimates by method');
    xtickangle(30);
    exportgraphics(fig, fullfile(outDir, 'fig_method_frequency_comparison.png'), 'Resolution', 180);
    close(fig);

    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 1100, 620]);
    bar(categorical(resultTable.method(validRows)), resultTable.energy_fraction(validRows));
    grid on;
    ylabel('Projection energy fraction on original signal');
    title('Energy fraction explained by each estimated frequency');
    xtickangle(30);
    exportgraphics(fig, fullfile(outDir, 'fig_method_energy_fraction.png'), 'Resolution', 180);
    close(fig);

    % ========== 修改从这里开始 ==========
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 1100, 680]);
    tiledlayout(2, 1, 'TileSpacing', 'compact');
    
    % 子图1: SSA去噪（保持原样，显示前2000个点）
    nexttile;
    plot(t(1:2000), x(1:2000), 'Color', [0.6,0.6,0.6], 'LineWidth', 0.7);
    hold on;
    plot(t(1:2000), xSsa(1:2000), 'Color', [0.1,0.45,0.72], 'LineWidth', 1.0);
    grid on;
    title('SSA denoising preview');
    xlabel('Time / s');
    ylabel('Amplitude');
    legend({'Original', 'SSA rank-2'}, 'Location', 'best');
    
    % 子图2: 随机共振输出（限制时间轴为0-10秒）
    nexttile;
    % 找出0-10秒内的索引
    idx10 = find(srTime <= 10, 1, 'last');
    if isempty(idx10)
        idx10 = min(2000, numel(srTime));
    end
    plot(srTime(1:idx10), srOutput(1:idx10), 'Color', [0.15,0.55,0.35], 'LineWidth', 0.9);
    grid on;
    title('Stochastic resonance output preview (0-10 s)');
    xlabel('Time / s');
    ylabel('SR output');
    xlim([0, 10]);
    exportgraphics(fig, fullfile(outDir, 'fig_denoising_enhancement_preview.png'), 'Resolution', 180);
    close(fig);
    % ========== 修改结束 ==========

    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 1100, 760]);
    tiledlayout(numel(acfSignals), 1, 'TileSpacing', 'compact');
    for k = 1:numel(acfSignals)
        tau = acfTimes{k};
        y = acfSignals{k};
        [~, fit, ~] = fitSinusoidAtFrequency(tau, y, fProj);
        show = tau <= min(12, max(tau));
        nexttile;
        plot(tau(show), y(show), 'Color', [0.56,0.56,0.56], 'LineWidth', 0.8);
        hold on;
        plot(tau(show), fit(show), 'Color', [0.75,0.12,0.18], 'LineWidth', 1.2);
        grid on;
        title(sprintf('Autocorrelation enhancement order %d', k));
        xlabel('Lag / s');
        ylabel('Normalized value');
        legend({'Enhanced series', '2 Hz sinusoidal projection'}, 'Location', 'best');
    end
    exportgraphics(fig, fullfile(outDir, 'fig_autocorrelation_enhancement.png'), 'Resolution', 180);
    close(fig);
end
function writeAnalysisReport(outDir, tbl, fs, duration, dfTheory, fFftBin, fProj, bestGain, srProj, srSnrDb)
    reportFile = fullfile(outDir, 'q1_method_comparison_report_auto.md');
    fid = fopen(reportFile, 'w', 'n', 'UTF-8');
    fprintf(fid, '# 第1问多方法频率检测结果分析\n\n');
    fprintf(fid, '采样频率为 %.6f Hz，观测时长为 %.6f s，全长理论频率分辨率为 %.12f Hz。\n\n', fs, duration, dfTheory);
    fprintf(fid, '普通全长 FFT 的最大频点为 %.12f Hz；本文以全长投影/GLRT 的 %.12f Hz 作为主参考结果。\n\n', fFftBin, fProj);
    fprintf(fid, '## 方法对比表\n\n');
    fprintf(fid, '| 方法 | 频率/Hz | 能量占比 | 与投影法差值/Hz | 评价 |\n');
    fprintf(fid, '|---|---:|---:|---:|---|\n');
    for i = 1:height(tbl)
        fprintf(fid, '| %s | %.12g | %.12g | %.12g | %s |\n', ...
            string(tbl.method{i}), tbl.frequency_Hz(i), tbl.energy_fraction(i), ...
            tbl.difference_from_projection_Hz(i), string(tbl.comment{i}));
    end
    fprintf(fid, '\n## 结论\n\n');
    fprintf(fid, '1. 全长投影/GLRT 与密集相关扫描本质一致，均利用 cos/sin 双参考信号进行相干积累，结果最稳定，适合作为第一问主方法。\n');
    fprintf(fid, '2. 插值 FFT 和零填充 FFT 能修正普通 FFT 的频点栅栏误差，但统计判决和参数恢复能力不如投影/GLRT 完整。\n');
    fprintf(fid, '3. AR 预白化可检验有色噪声影响。若其频率结果与主方法一致，可作为鲁棒性证据。\n');
    fprintf(fid, '4. SSA 去噪可作为数据驱动预处理，但可能改变幅值结构，适合作辅助展示，不建议替代主频率估计。\n');
    fprintf(fid, '5. 随机共振输出也能检测到接近 2 Hz 的频率，但增强优势有限。当前最佳输入增益 %.6g，输出投影占比 %.6g，局部 SNR %.6g dB，适合作非线性辅助验证。\n', bestGain, srProj, srSnrDb);
    fprintf(fid, '\n综合建议：第一问主线采用“全长投影搜索 + GLRT/相关检测”；插值 FFT、AR 预白化、SSA 和随机共振作为对比或鲁棒性分析。\n');
    fclose(fid);
end
