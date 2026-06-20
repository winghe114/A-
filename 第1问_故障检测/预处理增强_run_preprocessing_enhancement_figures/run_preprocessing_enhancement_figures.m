%% Preprocessing enhancement figures for Question 1
% Generate figures and data used by Section 5.1.1.

clear; clc; close all;

scriptDir = fileparts(mfilename('fullpath'));
q1Dir = fileparts(scriptDir);
rootDir = fileparts(q1Dir);
dataFile = fullfile(rootDir, 'data.xlsx');
outputDir = scriptDir;

fprintf('Generating preprocessing enhancement figures...\n');

[t, xRaw] = readSingleSourceData(dataFile);
t = t(:);
xRaw = xRaw(:);
y = xRaw - mean(xRaw, 'omitnan');
sigmaY = std(y, 0, 'omitnan');
u = y / sigmaY;
dt = median(diff(t));
fs = 1 / dt;
f0 = 2.0000016152054;

%% Stochastic resonance output
srParams = struct();
srParams.a = 1.0;
srParams.b = 1.0;
srParams.downsampleFactor = 2;
srParams.keepFraction = 0.10;
srParams.initialState = 0.05;
srParams.smoothingWindow = max(5, round(fs / srParams.downsampleFactor / f0 / 8));
bestGain = 0.2733;
[srTime, srOutput] = simulateSr(t, u, srParams, bestGain);

%% First-order autocorrelation enhancement
maxLagSeconds = 200;
[acfLag, acfSignal, acfRaw] = firstOrderAutocorrelation(u, fs, maxLagSeconds);

%% Save data
srTable = table(srTime, srOutput, 'VariableNames', {'time_s', 'sr_output'});
writetable(srTable, fullfile(outputDir, 'stochastic_resonance_output.csv'));

acfTable = table(acfLag, acfSignal, acfRaw, ...
    'VariableNames', {'lag_s', 'acf_standardized', 'acf_raw_biased'});
writetable(acfTable, fullfile(outputDir, 'autocorrelation_order1.csv'));

%% Figures
makeBistablePotentialFigure(outputDir, srParams.a, srParams.b);
makeSrComparisonFigure(outputDir, t, u, srTime, srOutput, fs, f0);
makeAcfComparisonFigure(outputDir, t, u, acfLag, acfSignal, fs, f0);

fprintf('Finished. Outputs saved in %s\n', outputDir);

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

function [lag, acfSignal, acfRaw] = firstOrderAutocorrelation(u, fs, maxLagSeconds)
    u = u(:) - mean(u(:), 'omitnan');
    n = numel(u);
    nfft = 2 ^ nextpow2(2 * n - 1);
    U = fft(u, nfft);
    r = real(ifft(U .* conj(U)));
    r = r(1:n) / n;
    r = r(2:end);
    keepCount = min(numel(r), max(10, floor(maxLagSeconds * fs)));
    acfRaw = r(1:keepCount);
    lag = (1:keepCount)' / fs;
    acfSignal = detrend(acfRaw, 'linear');
    acfSignal = normalizeColumn(acfSignal);
end

function y = normalizeColumn(x)
    y = x(:) - mean(x(:), 'omitnan');
    scale = std(y, 0, 'omitnan');
    if ~isfinite(scale) || scale <= eps
        scale = 1;
    end
    y = y / scale;
end

function makeBistablePotentialFigure(outDir, a, b)
    z = linspace(-2.2, 2.2, 1000)';
    U = -0.5 * a * z.^2 + 0.25 * b * z.^4;
    zStable = sqrt(a / b);
    barrier = a^2 / (4 * b);
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 980, 560]);
    plot(z, U, 'Color', [0.1,0.28,0.58], 'LineWidth', 2.0);
    hold on;
    scatter([-zStable, zStable], [-barrier, -barrier], 70, [0.75,0.12,0.18], 'filled');
    scatter(0, 0, 70, [0.12,0.55,0.32], 'filled');
    xline(-zStable, '--', 'Color', [0.75,0.12,0.18], 'LineWidth', 1.0);
    xline(zStable, '--', 'Color', [0.75,0.12,0.18], 'LineWidth', 1.0);
    yline(0, ':', 'Color', [0.25,0.25,0.25], 'LineWidth', 1.0);
    grid on;
    xlabel('System state z');
    ylabel('Potential U(z)');
    title('Bistable potential for stochastic resonance, a = b = 1');
    text(-zStable-0.55, -barrier-0.10, 'stable state z=-1', 'Color', [0.75,0.12,0.18]);
    text(zStable+0.08, -barrier-0.10, 'stable state z=1', 'Color', [0.75,0.12,0.18]);
    text(0.08, 0.08, 'barrier at z=0', 'Color', [0.12,0.55,0.32]);
    text(-2.05, 0.45, sprintf('Barrier height \\DeltaU = %.2f', barrier), ...
        'FontWeight', 'bold', 'Color', [0.1,0.1,0.1]);
    exportgraphics(fig, fullfile(outDir, 'fig_sr_bistable_potential.png'), 'Resolution', 220);
    close(fig);
end

function makeSrComparisonFigure(outDir, t, u, srTime, srOutput, fs, f0)
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 1180, 820]);
    tiledlayout(2, 1, 'TileSpacing', 'compact');
    tStart = srTime(1);
    tEnd = min(tStart + 12, srTime(end));
    showOriginal = t >= tStart & t <= tEnd;
    showSr = srTime >= tStart & srTime <= tEnd;
    srPlot = normalizeColumn(srOutput);
    nexttile;
    plot(t(showOriginal), u(showOriginal), 'Color', [0.58,0.58,0.58], 'LineWidth', 0.75);
    hold on;
    plot(srTime(showSr), srPlot(showSr), 'Color', [0.75,0.12,0.18], 'LineWidth', 1.15);
    grid on;
    xlabel('Time / s');
    ylabel('Normalized amplitude');
    title('Short-time waveform before and after stochastic resonance');
    legend({'Standardized original signal u(t)', 'Stochastic resonance output z(t)'}, 'Location', 'best');

    [fOrig, ampOrig] = amplitudeSpectrum(u, fs);
    fsSr = 1 / median(diff(srTime));
    [fSr, ampSr] = amplitudeSpectrum(srPlot, fsSr);
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
    title('Amplitude spectrum before and after stochastic resonance');
    legend({'Original spectrum', 'SR output spectrum', 'Estimated fault frequency'}, 'Location', 'best');
    exportgraphics(fig, fullfile(outDir, 'fig_sr_time_spectrum_comparison.png'), 'Resolution', 220);
    close(fig);
end

function makeAcfComparisonFigure(outDir, t, u, acfLag, acfSignal, fs, f0)
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 1180, 820]);
    tiledlayout(2, 1, 'TileSpacing', 'compact');
    showOriginal = t <= min(t(1) + 12, t(end));
    showAcf = acfLag <= min(12, max(acfLag));
    nexttile;
    plot(t(showOriginal), u(showOriginal), 'Color', [0.58,0.58,0.58], 'LineWidth', 0.75);
    hold on;
    plot(acfLag(showAcf), acfSignal(showAcf), 'Color', [0.1,0.45,0.72], 'LineWidth', 1.15);
    grid on;
    xlabel('Time or lag / s');
    ylabel('Normalized amplitude');
    title('Short-time waveform before and after first-order autocorrelation enhancement');
    legend({'Standardized original signal u(t)', 'Autocorrelation sequence R_u(\tau)'}, 'Location', 'best');

    [fOrig, ampOrig] = amplitudeSpectrum(u, fs);
    fsAcf = 1 / median(diff(acfLag));
    [fAcf, ampAcf] = amplitudeSpectrum(acfSignal, fsAcf);
    ampOrig = normalizeSpectrumInBand(fOrig, ampOrig, 0, 8);
    ampAcf = normalizeSpectrumInBand(fAcf, ampAcf, 0, 8);
    nexttile;
    plot(fOrig, ampOrig, 'Color', [0.58,0.58,0.58], 'LineWidth', 1.0);
    hold on;
    plot(fAcf, ampAcf, 'Color', [0.1,0.45,0.72], 'LineWidth', 1.2);
    xline(f0, '--', 'Color', [0.75,0.12,0.18], 'LineWidth', 1.2);
    xlim([0, 8]);
    grid on;
    xlabel('Frequency / Hz');
    ylabel('Normalized amplitude');
    title('Amplitude spectrum before and after autocorrelation enhancement');
    legend({'Original spectrum', 'Autocorrelation spectrum', 'Estimated fault frequency'}, 'Location', 'best');
    exportgraphics(fig, fullfile(outDir, 'fig_acf_time_spectrum_comparison.png'), 'Resolution', 220);
    close(fig);
end

function [freq, amp] = amplitudeSpectrum(x, fs)
    x = x(:) - mean(x(:), 'omitnan');
    n = numel(x);
    win = hannLocal(n);
    X = fft(x .* win);
    freq = (0:floor(n/2))' * fs / n;
    amp = abs(X(1:numel(freq))) / sum(win) * 2;
    if max(amp) <= eps
        amp(:) = 0;
    end
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

function w = hannLocal(n)
    k = (0:n-1)';
    w = 0.5 - 0.5 * cos(2*pi*k/max(n-1, 1));
end
