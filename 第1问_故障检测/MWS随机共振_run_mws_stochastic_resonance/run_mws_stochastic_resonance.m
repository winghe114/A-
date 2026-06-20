%% MWS stochastic resonance attempt for Question 1
% This script implements the MWS potential model described in Section 3.2.3
% and 3.2.4 of Shi Jiabei's thesis. It is used as an auxiliary enhancement
% experiment and does not replace the full-length sinusoidal projection model.

clear; clc; close all;

scriptDir = fileparts(mfilename('fullpath'));
q1Dir = fileparts(scriptDir);
rootDir = fileparts(q1Dir);
dataFile = fullfile(rootDir, 'data.xlsx');

fprintf('Running MWS stochastic resonance attempt...\n');

[t, xRaw] = readSingleSourceData(dataFile);
t = t(:);
xRaw = xRaw(:);
y = xRaw - mean(xRaw, 'omitnan');
u = y / std(y, 0, 'omitnan');
fs = 1 / median(diff(t));
f0 = 2.000001615205;

% Downsample and smooth the input for stochastic resonance integration. This
% follows the same practical treatment used by the previous SR branch.
baseParams = struct();
baseParams.a = 1.0;
baseParams.downsampleFactor = 2;
baseParams.keepFraction = 0.10;
baseParams.initialState = 0.05;
baseParams.stateLimit = 15;
baseParams.smoothingWindow = max(5, round(fs / baseParams.downsampleFactor / f0 / 8));

% Limited grid search. The MWS model has more parameters than the classical
% quartic bistable system, so this pass first uses cheap indicators.
bList = [0.8, 1.2, 1.8, 2.5, 3.5, 5.5];
v0List = [0.8, 1.5, 2.5, 4.0, 8.0, 16.0, 24.0];
rList = [0.25, 0.50, 0.80, 1.20, 1.60];
cList = [0.12, 0.20, 0.35, 0.60, 0.90];
kList = [0.04, 0.08, 0.14, 0.22, 0.32, 0.48, 0.70];

rows = {};
outputs = {};
times = {};
rowIdx = 0;

for ib = 1:numel(bList)
    for iv = 1:numel(v0List)
        for ir = 1:numel(rList)
            for ic = 1:numel(cList)
                for ik = 1:numel(kList)
                    p = baseParams;
                    p.b = bList(ib);
                    p.V0 = v0List(iv);
                    p.R = rList(ir);
                    p.C = cList(ic);
                    p.k = kList(ik);

                    [tm, xm, stableFlag] = simulateMwsSr(t, u, p);
                    if ~stableFlag || numel(xm) < 200
                        continue;
                    end

                    projFrac = coherentProjectionFraction(tm, xm, f0);
                    [localSnrDb, peakAmp] = localFrequencySnr(tm, xm, f0);
                    roughFreq = localPeakFrequency(tm, xm, f0, 0.35);

                    rowIdx = rowIdx + 1;
                    rows(rowIdx, :) = {p.a, p.b, p.V0, p.R, p.C, p.k, ...
                        projFrac, localSnrDb, peakAmp, roughFreq, NaN, NaN, NaN, NaN, NaN, NaN};
                    times{rowIdx, 1} = tm;
                    outputs{rowIdx, 1} = xm;
                end
            end
        end
    end
end

if isempty(rows)
    error('No stable MWS stochastic resonance output was obtained. Try smaller gains or wider state limits.');
end

resultTable = cell2table(rows, 'VariableNames', {'a', 'b', 'V0', 'R', 'C', 'k', ...
    'mws_projection_fraction', 'mws_local_snr_dB', 'mws_peak_amplitude', ...
    'mws_rough_frequency_Hz', 'mws_frequency_Hz', 'frequency_error_Hz', ...
    'abs_frequency_error_Hz', 'original_projection_energy', ...
    'original_energy_fraction', 'combined_score'});

% Score cheap indicators first, then refine the most promising parameter sets.
zProj = robustZ(resultTable.mws_projection_fraction);
zSnr = robustZ(resultTable.mws_local_snr_dB);
roughErr = abs(resultTable.mws_rough_frequency_Hz - f0);
zErr = -robustZ(roughErr);
quickScore = 0.42 * zProj + 0.43 * zSnr + 0.15 * zErr;
[~, order] = sort(quickScore, 'descend');
refineCount = min(24, numel(order));
refineIdx = order(1:refineCount);

totalEnergyOriginal = sum(y .^ 2);
for n = 1:numel(refineIdx)
    idx = refineIdx(n);
    tm = times{idx};
    xm = outputs{idx};
    [fHat, ~, ~, ~] = projectionSearchLight(tm, xm, f0, 0.35, 901);
    [~, fitOrig, ~] = fitSinusoidAtFrequency(t, y, fHat);
    origProjectionEnergy = sum(fitOrig .^ 2);
    origEnergyFraction = origProjectionEnergy / max(totalEnergyOriginal, eps);

    resultTable.mws_frequency_Hz(idx) = fHat;
    resultTable.frequency_error_Hz(idx) = fHat - f0;
    resultTable.abs_frequency_error_Hz(idx) = abs(fHat - f0);
    resultTable.original_projection_energy(idx) = origProjectionEnergy;
    resultTable.original_energy_fraction(idx) = origEnergyFraction;
end

refinedMask = isfinite(resultTable.mws_frequency_Hz);
if any(refinedMask)
    errScore = -robustZ(resultTable.abs_frequency_error_Hz(refinedMask));
    origScore = robustZ(resultTable.original_energy_fraction(refinedMask));
    projScore = robustZ(resultTable.mws_projection_fraction(refinedMask));
    snrScore = robustZ(resultTable.mws_local_snr_dB(refinedMask));
    refinedScore = 0.30 * projScore + 0.30 * snrScore + 0.25 * origScore + 0.15 * errScore;
    resultTable.combined_score(refinedMask) = refinedScore;
else
    resultTable.combined_score = quickScore;
end

[~, bestIdx] = max(resultTable.combined_score);
bestParams = table2struct(resultTable(bestIdx, :));
bestTime = times{bestIdx};
bestOutput = outputs{bestIdx};

% Classical quartic SR baseline using previously selected gain.
classicParams = struct();
classicParams.a = 1.0;
classicParams.b = 1.0;
classicParams.downsampleFactor = baseParams.downsampleFactor;
classicParams.keepFraction = baseParams.keepFraction;
classicParams.initialState = baseParams.initialState;
classicParams.stateLimit = baseParams.stateLimit;
classicParams.smoothingWindow = baseParams.smoothingWindow;
classicParams.k = 0.2733;
[classicTime, classicOutput] = simulateClassicSr(t, u, classicParams);
classicProj = coherentProjectionFraction(classicTime, classicOutput, f0);
[classicSnr, ~] = localFrequencySnr(classicTime, classicOutput, f0);
[classicFreq, ~, ~, ~] = projectionSearchLight(classicTime, classicOutput, f0, 0.35, 901);

% Save outputs.
writetable(resultTable, fullfile(scriptDir, 'mws_sr_parameter_scan_results.csv'));
bestOutputTable = table(bestTime, bestOutput, 'VariableNames', {'time_s', 'mws_output'});
writetable(bestOutputTable, fullfile(scriptDir, 'mws_sr_best_output.csv'));

makeMwsPotentialFigure(scriptDir, bestParams);
makeMwsComparisonFigure(scriptDir, t, u, fs, bestTime, bestOutput, f0, bestParams);
makeParameterScanFigure(scriptDir, resultTable, bestIdx, f0);
makeMwsVsClassicFigure(scriptDir, bestTime, bestOutput, classicTime, classicOutput, f0);
writeSummary(scriptDir, bestParams, f0, classicFreq, classicProj, classicSnr, resultTable);

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

function [tout, xout, stableFlag] = simulateMwsSr(t, xNorm, p)
    ds = max(1, round(p.downsampleFactor));
    td = t(1:ds:end);
    u = xNorm(1:ds:end);
    if p.smoothingWindow > 1
        u = movmean(u, p.smoothingWindow);
    end
    dt = median(diff(td));
    x = p.initialState;
    keepStart = max(1, floor(numel(td) * p.keepFraction));
    xout = zeros(numel(td)-keepStart, 1);
    tout = zeros(numel(td)-keepStart, 1);
    outIdx = 0;
    stableFlag = true;
    for n = 1:numel(td)-1
        inputValue = p.k * u(n);
        k1 = mwsDeriv(x, inputValue, p);
        k2 = mwsDeriv(x + 0.5*dt*k1, inputValue, p);
        k3 = mwsDeriv(x + 0.5*dt*k2, inputValue, p);
        k4 = mwsDeriv(x + dt*k3, inputValue, p);
        x = x + dt * (k1 + 2*k2 + 2*k3 + k4) / 6;
        if ~isfinite(x)
            stableFlag = false;
            break;
        end
        if abs(x) > p.stateLimit
            x = sign(x) * p.stateLimit;
        end
        if n >= keepStart
            outIdx = outIdx + 1;
            tout(outIdx) = td(n);
            xout(outIdx) = x;
        end
    end
    tout = tout(1:outIdx);
    xout = xout(1:outIdx);
    xout = xout - mean(xout, 'omitnan');
    if std(xout, 0, 'omitnan') <= eps
        stableFlag = false;
    end
end

function dx = mwsDeriv(x, inputValue, p)
    expArg = (abs(x) - p.R) / p.C;
    expArg = min(max(expArg, -60), 60);
    q = exp(expArg);
    wsForce = (p.V0 / p.C) * sign(x) * q / (1 + q)^2;
    dx = -p.b * x - wsForce + inputValue;
end

function [tout, yout] = simulateClassicSr(t, xNorm, p)
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
    for n = 1:numel(td)-1
        inputValue = p.k * u(n);
        k1 = classicDeriv(y, inputValue, p);
        k2 = classicDeriv(y + 0.5*dt*k1, inputValue, p);
        k3 = classicDeriv(y + 0.5*dt*k2, inputValue, p);
        k4 = classicDeriv(y + dt*k3, inputValue, p);
        y = y + dt * (k1 + 2*k2 + 2*k3 + k4) / 6;
        if abs(y) > p.stateLimit
            y = sign(y) * p.stateLimit;
        end
        if n >= keepStart
            outIdx = outIdx + 1;
            tout(outIdx) = td(n);
            yout(outIdx) = y;
        end
    end
    tout = tout(1:outIdx);
    yout = yout(1:outIdx);
    yout = yout - mean(yout, 'omitnan');
end

function dy = classicDeriv(y, inputValue, p)
    dy = p.a * y - p.b * y^3 + inputValue;
end

function frac = coherentProjectionFraction(t, y, f)
    y = y(:) - mean(y(:), 'omitnan');
    t = t(:);
    H = [cos(2*pi*f*t), sin(2*pi*f*t)];
    theta = H \ y;
    fit = H * theta;
    frac = sum(fit .^ 2) / max(sum(y .^ 2), eps);
end

function [snrDb, peakAmp] = localFrequencySnr(t, y, f0)
    y = y(:) - mean(y(:), 'omitnan');
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

function fPeak = localPeakFrequency(t, y, f0, halfWidth)
    fs = 1 / median(diff(t));
    [freq, amp] = amplitudeSpectrum(y, fs);
    band = freq >= f0 - halfWidth & freq <= f0 + halfWidth;
    fb = freq(band);
    ab = amp(band);
    [~, idx] = max(ab);
    fPeak = fb(idx);
end

function [theta, fit, residual] = fitSinusoidAtFrequency(t, x, f)
    t = t(:);
    x = x(:) - mean(x(:), 'omitnan');
    H = [cos(2*pi*f*t), sin(2*pi*f*t)];
    theta = H \ x;
    fit = H * theta;
    residual = x - fit;
end

function [fHat, coarseFreq, projEnergy, residualEnergy] = projectionSearchLight(t, x, center, halfWidth, gridN)
    t = t(:);
    x = x(:) - mean(x(:), 'omitnan');
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
    w = hannLocal(n);
    X = fft(x .* w);
    freq = (0:floor(n/2))' * fs / n;
    amp = abs(X(1:numel(freq))) / sum(w) * 2;
end

function w = hannLocal(n)
    k = (0:n-1)';
    w = 0.5 - 0.5 * cos(2*pi*k/max(n-1, 1));
end

function z = robustZ(x)
    x = x(:);
    finiteMask = isfinite(x);
    z = zeros(size(x));
    if nnz(finiteMask) < 2
        return;
    end
    mu = median(x(finiteMask), 'omitnan');
    sig = 1.4826 * median(abs(x(finiteMask) - mu), 'omitnan');
    if ~isfinite(sig) || sig <= eps
        sig = std(x(finiteMask), 0, 'omitnan');
    end
    if ~isfinite(sig) || sig <= eps
        sig = 1;
    end
    z(finiteMask) = (x(finiteMask) - mu) / sig;
end

function y = normalizeColumn(x)
    y = x(:) - mean(x(:), 'omitnan');
    s = std(y, 0, 'omitnan');
    if ~isfinite(s) || s <= eps
        s = 1;
    end
    y = y / s;
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

function U = mwsPotential(x, p)
    expArg = (abs(x) - p.R) / p.C;
    expArg = min(max(expArg, -60), 60);
    U = -p.a + 0.5 * p.b * x.^2 + p.V0 ./ (1 + exp(expArg));
end

function makeMwsPotentialFigure(outDir, p)
    x = linspace(-3, 3, 1600)';
    U = mwsPotential(x, p);
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 980, 560]);
    plot(x, U, 'Color', [0.12,0.34,0.62], 'LineWidth', 2.0);
    grid on;
    xlabel('System state x');
    ylabel('Potential U(x)');
    title(sprintf('MWS potential, b=%.3g, V_0=%.3g, R=%.3g, C=%.3g', p.b, p.V0, p.R, p.C));
    exportgraphics(fig, fullfile(outDir, 'fig_mws_potential.png'), 'Resolution', 220);
    close(fig);
end

function makeMwsComparisonFigure(outDir, t, u, fs, tm, xm, f0, p)
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 1180, 820]);
    tiledlayout(2, 1, 'TileSpacing', 'compact');
    tStart = tm(1);
    tEnd = min(tStart + 12, tm(end));
    showOriginal = t >= tStart & t <= tEnd;
    showMws = tm >= tStart & tm <= tEnd;
    xmPlot = normalizeColumn(xm);

    nexttile;
    plot(t(showOriginal), u(showOriginal), 'Color', [0.58,0.58,0.58], 'LineWidth', 0.75);
    hold on;
    plot(tm(showMws), xmPlot(showMws), 'Color', [0.78,0.16,0.18], 'LineWidth', 1.15);
    grid on;
    xlabel('Time / s');
    ylabel('Normalized amplitude');
    title(sprintf('MWS SR output, b=%.3g, V_0=%.3g, R=%.3g, C=%.3g, k=%.3g', ...
        p.b, p.V0, p.R, p.C, p.k));
    legend({'Standardized original signal u(t)', 'MWS SR output'}, 'Location', 'best');

    [fOrig, ampOrig] = amplitudeSpectrum(u, fs);
    fsMws = 1 / median(diff(tm));
    [fMws, ampMws] = amplitudeSpectrum(xmPlot, fsMws);
    ampOrig = normalizeSpectrumInBand(fOrig, ampOrig, 0, 8);
    ampMws = normalizeSpectrumInBand(fMws, ampMws, 0, 8);
    nexttile;
    plot(fOrig, ampOrig, 'Color', [0.58,0.58,0.58], 'LineWidth', 1.0);
    hold on;
    plot(fMws, ampMws, 'Color', [0.78,0.16,0.18], 'LineWidth', 1.2);
    xline(f0, '--', 'Color', [0.10,0.30,0.58], 'LineWidth', 1.2);
    xlim([0, 8]);
    grid on;
    xlabel('Frequency / Hz');
    ylabel('Normalized amplitude');
    title('Spectrum before and after MWS stochastic resonance');
    legend({'Original spectrum', 'MWS SR spectrum', 'Full-length projection frequency'}, 'Location', 'best');
    exportgraphics(fig, fullfile(outDir, 'fig_mws_time_spectrum_comparison.png'), 'Resolution', 220);
    close(fig);
end

function makeParameterScanFigure(outDir, tbl, bestIdx, f0)
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 1180, 900]);
    tiledlayout(3, 1, 'TileSpacing', 'compact');
    x = 1:height(tbl);
    nexttile;
    plot(x, tbl.mws_local_snr_dB, '.', 'Color', [0.75,0.12,0.18], 'MarkerSize', 8);
    hold on;
    scatter(bestIdx, tbl.mws_local_snr_dB(bestIdx), 70, [0.08,0.32,0.62], 'filled');
    grid on;
    xlabel('Parameter combination index');
    ylabel('Local SNR / dB');
    title('MWS parameter scan: local SNR around 2 Hz');

    nexttile;
    plot(x, tbl.mws_projection_fraction, '.', 'Color', [0.10,0.52,0.36], 'MarkerSize', 8);
    hold on;
    scatter(bestIdx, tbl.mws_projection_fraction(bestIdx), 70, [0.08,0.32,0.62], 'filled');
    grid on;
    xlabel('Parameter combination index');
    ylabel('Projection fraction');
    title('MWS parameter scan: output projection fraction');

    nexttile;
    freqForPlot = tbl.mws_rough_frequency_Hz;
    refined = isfinite(tbl.mws_frequency_Hz);
    freqForPlot(refined) = tbl.mws_frequency_Hz(refined);
    plot(x, abs(freqForPlot - f0), '.', 'Color', [0.30,0.30,0.30], 'MarkerSize', 8);
    hold on;
    scatter(bestIdx, abs(freqForPlot(bestIdx) - f0), 70, [0.08,0.32,0.62], 'filled');
    grid on;
    xlabel('Parameter combination index');
    ylabel('|frequency error| / Hz');
    title('MWS parameter scan: frequency deviation from full-length projection');

    exportgraphics(fig, fullfile(outDir, 'fig_mws_parameter_scan.png'), 'Resolution', 220);
    close(fig);
end

function makeMwsVsClassicFigure(outDir, tm, xm, tc, xc, f0)
    xmPlot = normalizeColumn(xm);
    xcPlot = normalizeColumn(xc);
    fsMws = 1 / median(diff(tm));
    fsClassic = 1 / median(diff(tc));
    [fMws, ampMws] = amplitudeSpectrum(xmPlot, fsMws);
    [fClassic, ampClassic] = amplitudeSpectrum(xcPlot, fsClassic);
    ampMws = normalizeSpectrumInBand(fMws, ampMws, 0, 8);
    ampClassic = normalizeSpectrumInBand(fClassic, ampClassic, 0, 8);

    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 1180, 620]);
    plot(fMws, ampMws, 'Color', [0.78,0.16,0.18], 'LineWidth', 1.25);
    hold on;
    plot(fClassic, ampClassic, 'Color', [0.12,0.45,0.72], 'LineWidth', 1.25);
    xline(f0, '--', 'Color', [0.10,0.10,0.10], 'LineWidth', 1.2);
    xlim([0, 8]);
    grid on;
    xlabel('Frequency / Hz');
    ylabel('Normalized amplitude');
    title('MWS stochastic resonance vs classical quartic stochastic resonance');
    legend({'MWS SR output', 'Classical SR output', 'Full-length projection frequency'}, 'Location', 'best');
    exportgraphics(fig, fullfile(outDir, 'fig_mws_vs_classic_sr.png'), 'Resolution', 220);
    close(fig);
end

function writeSummary(outDir, best, f0, classicFreq, classicProj, classicSnr, tbl)
    summaryFile = fullfile(outDir, 'mws_sr_summary.txt');
    fid = fopen(summaryFile, 'w');
    cleaner = onCleanup(@() fclose(fid));
    fprintf(fid, 'MWS 型随机共振尝试结果\n');
    fprintf(fid, '========================\n\n');
    fprintf(fid, '主模型全长正弦投影频率: %.12f Hz\n\n', f0);
    fprintf(fid, 'MWS 势函数: U(x) = -a + b*x^2/2 + V0/(1+exp((|x|-R)/C))\n');
    fprintf(fid, 'MWS 系统方程: dx/dt = -b*x - (V0/C)*sgn(x)*q/(1+q)^2 + k*u(t), q=exp((|x|-R)/C)\n\n');
    fprintf(fid, '扫描参数组合数: %d\n', height(tbl));
    fprintf(fid, '最优参数: a=%.6g, b=%.6g, V0=%.6g, R=%.6g, C=%.6g, k=%.6g\n', ...
        best.a, best.b, best.V0, best.R, best.C, best.k);
    fprintf(fid, 'MWS 输出投影能量占比: %.12g\n', best.mws_projection_fraction);
    fprintf(fid, 'MWS 输出局部信噪比: %.6f dB\n', best.mws_local_snr_dB);
    fprintf(fid, 'MWS 粗频率估计: %.12f Hz\n', best.mws_rough_frequency_Hz);
    if isfinite(best.mws_frequency_Hz)
        fprintf(fid, 'MWS 精细投影频率估计: %.12f Hz\n', best.mws_frequency_Hz);
        fprintf(fid, '与主模型频率差: %.12g Hz\n', best.frequency_error_Hz);
        fprintf(fid, '与主模型绝对频率差: %.12g Hz\n', best.abs_frequency_error_Hz);
        fprintf(fid, '将 MWS 频率回代原始信号后的投影能量占比: %.12g\n', best.original_energy_fraction);
    end
    fprintf(fid, '\n经典四次双稳态 SR 对比:\n');
    fprintf(fid, '经典 SR 频率估计: %.12f Hz\n', classicFreq);
    fprintf(fid, '经典 SR 输出投影能量占比: %.12g\n', classicProj);
    fprintf(fid, '经典 SR 输出局部信噪比: %.6f dB\n\n', classicSnr);
    fprintf(fid, '结论:\n');
    fprintf(fid, '本文尝试了更复杂的 MWS 型随机共振势函数。若 MWS 输出在 2 Hz 附近有增强，但频率精度或回代原始信号解释能力仍不优于全长正弦投影，则 MWS 仅作为辅助增强实验，不参与最终频率估计。\n');
end
